# Authors: Leo McHugh, Immunexpress Inc
#          Kevin Snyder, FDA
# Code for generating simulation datasets
# Code for introducing noise into loaded datasets

library(shiny) # for GUI
library(caret) # for confusion matrices
library(pROC) # for ROC
library(e1071)
library(shinycssloaders)

#################################################################
######################## FUNCTIONS ##############################
#################################################################
# bounds a list of values between upper and lower limits
# used to truncate values (ie probabilities) between zero and one. 
bound = function(vals,lower,upper) {  
  output = vals;
  for (i in 1:length(output)) {
    if (output[i] < lower) {output[i]=lower}
    if (output[i] > upper) {output[i]=upper}
  }
  return(output)
}

# given a list of values and a threshold, returns a list of calls based on the thresh
# if a value is above the threshold, a value of 1 is returned, else, a value of zero
# corresponding to positive and negative test calls for a condition
thresh.calls = function(test.values,thresh) {
  binary.calls = c()
  for (val in test.values) {
    if (is.na(val)) {binary.calls=c(binary.calls,NA); next}
    if (val < thresh) {binary.calls=c(binary.calls,0); next} # zero for negative
    binary.calls=c(binary.calls,1) # else one for positive
  }
  return(binary.calls)
}

# return the median and empirical 95% CI values for a vector of values
med.CI = function(vals) {
  med = median(vals)
  CI.lo = sort(vals)[floor(length(vals)*0.05)]
  CI.hi = sort(vals)[ceiling(length(vals)*0.95)]
  output = c(med,CI.lo,CI.hi)
  names(output)=c("median",'CI.lo','CI.hi')
  return(output)
}

# Takes in a list of values and returns the median and 95% CI - same as above function? Remove?
condense.CI = function(med.lo.hi) {
  return(paste(sprintf("%0.3f",med.lo.hi[1]),' (',
               sprintf("%0.3f",med.lo.hi[2]),'-',
               sprintf("%0.3f",med.lo.hi[3]),')',sep=''))
}

# given an overall missclassification rate
# and a number of positive and negative ground truths
# and a different weight of probability of misclassification for positive and negative (arbitrary)
# solve for the scaling factor that gives positive and negative misclassification rates to result in a
# fixed overall misclassification rate. 
# mc = misclassification ; wt = weight
scaling.factor = function(overall.mc.rate,n.pos,pos.mc.wt,n.neg,neg.mc.wt) {
  # n.pos=200; n.neg=100; pos.mc.wt=2; neg.mc.wt=1; overall.mc.rate=0.10
  scale.factor = ((n.pos+n.neg)*overall.mc.rate) / (n.pos*pos.mc.wt + n.neg*neg.mc.wt)
  # out = c(pos.mc.wt,neg.mc.wt)*scale.factor; names(out)=c('neg.misclass.rate','pos.misclass.rate') !!! THIS IS MIXED UP !!!
  out = c(pos.mc.wt,neg.mc.wt)*scale.factor; names(out)=c('pos.misclass.rate','neg.misclass.rate')
  return(out)
  # pos.mc.wt*n.pos*scale.factor+neg.mc.wt*n.neg*scale.factor # debugging
}

# given a list of numbers, left pad them to accommodate the largest
# ie from 1 to 100: 001, 002 ... 100
my.pad = function(values) {
  pad.digits = nchar(as.character(max(values)))
  input.cmd = paste('%0',pad.digits,'.0f',sep='')
  return(sapply(values,function(val){return(sprintf(input.cmd, val))}))
}

# helper function to format the performance specs object (custom) into a string format
little.perf.table = function(given.ROC,given.cm) {
  return(
    paste(
      paste('AUC:',round(given.ROC$auc,digits=3)),
      paste('Sens:',sprintf("%0.3f",given.cm$byClass['Sensitivity'])),
      paste('Spec:',sprintf("%0.3f",given.cm$byClass['Specificity'])),
      paste('PPV :',sprintf("%0.3f",given.cm$byClass['Pos Pred Value'])), 
      paste('NPV :',sprintf("%0.3f",given.cm$byClass['Neg Pred Value'])),
      
      sep='\n'
      
    )
  )
}

# custom utility function for plotting 
plot.data.point = function(x,y,CI.lo,CI.hi,h.offset,xlim,ylim,color) {
  #
  #x=pcs[i]; y=point['median'];CI.lo = point['CI.lo']; CI.hi=point['CI.hi']; color=metric.cols[i]
  #
  par(new=T)
  plot(x+h.offset,y,axes=F,xlim=xlim,ylim=ylim,xlab='',ylab='',pch=16,col=color)
  segments(x+h.offset,CI.lo,x+h.offset,CI.hi,col=color,lwd=2)
  segments(x+h.offset-0.3,CI.hi,x+h.offset+0.3,CI.hi,col=color,lwd=2)
  segments(x+h.offset-0.3,CI.lo,x+h.offset+0.3,CI.lo,col=color,lwd=2)
}

#  ---- function introducing error ---- 
# ref.probs contains reference probabilities, the seed is for random consistency. 
create.reference = function(ref.probs,seed) {
  dirty.ref = c() # not clean - contains errors randomly sampled from the ref.probs. 
  for (i in 1:length(ref.probs)) {
    dirty.ref = c(dirty.ref,sample(c(1,0),1,prob=c(ref.probs[i],(1-ref.probs[i]))))
  }
  return(dirty.ref)
}

# produces as table formatted as text
table.performance.results = function(performance.results) {
  text.results = apply(performance.results,2,function(this.col){return(condense.CI(med.CI(this.col)))})
  cats = c('True','FV','MR','Adjusted')
  cats.found = c()
  metrics = c('AUC','Sensitivity','Specificity','Pos.Pred.Value','Neg.Pred.Value')
  metrics.found = c()
  for (m in metrics) {if (length(grep(m,names(text.results)))!=0) {metrics.found=c(metrics.found,m)}}
  out.table=c()
  for (cat in cats) {
    this.cat = grep(cat,names(text.results))
    if (length(this.cat)==0) {next}
    cats.found = c(cats.found,cat)
    out.table=cbind(out.table,text.results[this.cat])
  }; colnames(out.table)=cats.found; rownames(out.table)=metrics.found
  return(out.table)
}

server <- function(input, output,session) {
  
  #########################################################
  ## SIMULATIONS ON THE EFFECT OF INCREASING UNCERTAINTY ##
  #########################################################
  
  # ===================== RUN SIMULATION ============================
  
  runSimulation <- eventReactive(input$go,{
    withProgress(message = '',value=0, expr={
      
      # Update Parameter Values
      n <- input$n
      p <- input$p
      sens <- input$sens
      spec <- input$spec
      pos.mc.wt <- input$FP
      neg.mc.wt <- input$FN
      reps <- 100
      
      # Create Simulated Data Set
      n.pos <- round(n*p,digits=0)
      n.neg <- round(n*(1-p),digits=0)
      gtPos <- rep(1,n.pos)
      gtNeg <- rep(0,n.neg)
      ground.truth = c(gtPos,gtNeg) # ground truth
      gtPosTestResult <- c(rep(1,round(n.pos*sens)),rep(0,n.pos-round(n.pos*sens)))
      gtNegTestResult <- c(rep(0,round(n.neg*spec)),rep(1,n.neg-round(n.neg*spec)))
      Dx.call <- c(gtPosTestResult,gtNegTestResult)
      
      # Simulate Comparator Data Sets
      mc.rates =  seq(from=0,to=0.5,by=0.05) # misclassification rates
      results=list()
      mc.wts.rates <- cbind(NULL,NULL)
      for (misclass.rate in mc.rates) {
        this.rate = c()
        mc.wts = scaling.factor(misclass.rate,n.pos,pos.mc.wt,n.neg,neg.mc.wt) 
        mc.wts.rates <- rbind(mc.wts.rates,mc.wts)
        for (my.seed in seq(reps)) {
          noisy.probs = ground.truth; noisy.probs[which(noisy.probs==1)]=(1-mc.wts['pos.misclass.rate']); noisy.probs[which(noisy.probs==0)]=mc.wts['neg.misclass.rate']
          reference = create.reference(noisy.probs,my.seed)
          ROC = pROC::roc(predictor=Dx.call,response=reference,ci=T)
          cm = confusionMatrix(factor(Dx.call),factor(reference),positive = '1')
          this.row = c(misclass.rate,ROC$auc,cm$byClass['Sensitivity'],cm$byClass['Specificity'],cm$byClass['Pos Pred Value'],cm$byClass['Neg Pred Value'])
          names(this.row)[1]='misclass.rate'; names(this.row)[2]='AUC'
          this.rate = rbind(this.rate,this.row)
        }
        results = append(results,list(this.rate))
        setProgress(message='Simulation in Progress...',value=misclass.rate/tail(mc.rates,1))
      }
      names(results)=as.character(mc.rates)
      outList <- list()
      outList$results <- results
      outList$mc.wts.rates <- mc.wts.rates
      return(outList)
    })
  })
  
  
  # ======================== PLOT OUT ===============================
  
  output$plot <- renderPlot({
    outList <- runSimulation()
    results <- outList$results
    mc.wts.rates <- outList$mc.wts.rates
    row.names(mc.wts.rates) <- seq(nrow(mc.wts.rates))
    mc.wts.rates <- as.data.frame(mc.wts.rates)
    
    # set up the background
    par(mar=c(7,5,1,10))
    par(xpd=T)
    xlim=c(-1,51); ylim=c(0,1)#perf.ylim
    par(new=F)
    plot(NA,xlim=xlim,ylim=ylim,axes=F,xlab='',ylab='',cex.lab=1.4)
    axis(1,at=c(0,5,10,15,20,25,30,35,40,45,50),labels=round(mc.wts.rates$pos.misclass.rate*100,digits=1),cex.axis=1.5)
    axis(1,at=c(0,5,10,15,20,25,30,35,40,45,50),labels=round(mc.wts.rates$neg.misclass.rate*100,digits=1),line=3,cex.axis=1.5)
    axis(2,at=seq(from=0.0,to=1,by=0.1),cex.axis=1.5,las=2)
    mtext('Comparator Misclassification rate (%)',side=1,line=6,cex=1.5)
    mtext('FP rate:',1,line=0.5,at=-3.5,cex=1.5)
    mtext('FN rate:',1,line=3.5,at=-3.5,cex=1.5)
    for (lvl in seq(from=0.0,to=1,by=0.1)) {abline(h=lvl,lty=2,xpd=F)} # horiz grid lines
    pcs = c(0,5,10,15,20,25,30,35,40,45,50)
    metrics =  c("AUC","Sensitivity","Specificity","Pos Pred Value","Neg Pred Value") #
    offsets = c(-0.8,-0.4,0,0.4,0.8); names(offsets)=c("AUC","Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")
    metric.cols = c('red','blue','black','darkgreen','darkmagenta'); names(metric.cols)=c("AUC","Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")
    
    # populate the plot
    for (metric in metrics) {
      prev.point=NA
      for (i in 1:length(results)) {
        vals = results[[i]][,metric]
        point = med.CI(vals)
        offset = offsets[metric]
        plot.data.point(pcs[i],point['median'],point['CI.lo'],point['CI.hi'],offset,xlim,ylim,metric.cols[metric])
        this.point = c(pcs[i]+offset,point['median'])
        if (length(prev.point)==1) {prev.point=this.point; next} # if NA
        segments(this.point[1],this.point[2],prev.point[1],prev.point[2],col=metric.cols[metric],lwd=2); prev.point=this.point;
      }
    }
    legend(55,1,legend=c("AUC ","PPA ","NPA ","PPV ","NPV "),col=metric.cols,pch=16,lwd=1,cex=1.5)
  })
  
}

############################### Define GUI for Application #################################

ui <- fluidPage(
  
  titlePanel("The Effect of Uncertainty in Patient Classification on Diagnostic Performance Estimation"),
  
  sidebarLayout(
    
    sidebarPanel(
      numericInput('n','Number of Patients:',200,step=1),
      sliderInput('p','Disease Prevalence:',0,1,0.5,step=0.01),
      sliderInput('sens','True Test Sensitivity:',0,1,0.9,step=0.01),
      sliderInput('spec','True Test Specificity:',0,1,0.9,step=0.01),
      h4('Comparator False Positive (FP) to False Negative (FN) Ratio*'),
      splitLayout(
        shiny::div(style="text-align:center",numericInput('FP','FP',1,step=1,width='60px')),
        HTML('<div style="line-height:200%;"><br></div><big>:</big>'),
        shiny::div(style="text-align:center",numericInput('FN','FN',1,step=1,width='60px')),
        cellWidths=c('60px','5px','60px')
      ),
      HTML('* For example, a ratio of 1:1 means that misclassifications are equally likely for positive and negative samples whereas
         a ratio of 2:1 means that true negative samples are twice as likely to be misclassified as positive (FP)
         than true positive samples are likely to be misclassified as negative (FN).'),br(),br(),
      actionButton('go','Perform Simulation',icon=icon('arrow-circle-right')),
      br(),br(),
      actionButton('source','View Source Code',icon=icon('code'),
                   onclick='window.open("https://github.com/ksny/Imperfect-Gold-Standard/blob/master/MakingNoiseApp.R","_blank")')
    ),
    
    mainPanel(
      withSpinner(plotOutput('plot',height = 600),type=5),#,color.background='white'),
      conditionalPanel(condition = 'input.go > 0',hr(),
                       HTML('<P ALIGN=RIGHT>AUC = Area Under Curve<br>PPA = Positive Percent Agreement<br>NPA = Negative Percent Agreement<br>
                            PPV = Positive Predictive Value<br>NPV = Negative Predictive Value</P>')
      )
    )
  )
)

############################################################################################

# Run Shiny App
shinyApp(ui = ui, server = server)



