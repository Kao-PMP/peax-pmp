library(shiny)
library(parallel)
library(miscTools)
#library(multiMiR)
source('df2html.R')

oldselmir=0
oldselmir1=0
oldselmir2=0
oldsearch=0
oldsearchall=0
oldloaddemo=0
demoloads=0
numtests=0
curfilt=0

shinyServer(function(input, output) {
  # output$tablespacing <- renderText({
  #   "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n."
  # })
  # output$tablespacing2 <- renderText({
  #   "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n."
  # })
  
  # testcor2<-reactive({
  #   return(testcor(2))
  # })
  # 
  # testcor1<-reactive({
  #   return(testcor(1))
  # })
  
  datasetInput<-reactive({
    inFile <- input$file1
    
    if (!is.null(inFile)) {
      #message("READING DATASET")
      #message(inFile$datapath)
      x<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    }
    else if (input$loaddemo > 0) {
      #message("DEFAULT DATASET")
      x<-read.csv("./data/peax_demo_pheno.csv", header=input$header, sep=input$sep, quote=input$quote)
      #x<-read.csv("./data/combined-5.csv", header=input$header, sep=input$sep, quote=input$quote)
      #x<-read.csv("./data/pilot_4.csv", header=input$header, sep=input$sep, quote=input$quote)
    }
    else
    {
      return(NULL) }
    return(x)})
  
  # mirInput<-reactive({
  #   #input$search
  #   #message("READING MIRINPUT")
  #   mirFile <- input$mirfile
  #   
  #   if (!is.null(mirFile))
  #     y<-read.csv(mirFile$datapath,check.names=FALSE,header=input$header, sep=input$sep, quote=input$quote)
  #   else if (input$loaddemo > 0) {
  #     #message("DEFAULT DATASET")
  #     y<-read.csv("./data/peax_demo_mrna2.csv", header=TRUE, sep=",", quote="\"")
  #   }
  #   else
  #     return(NULL)
  #   
  #   if (nrow(y)>50000) {
  #     y<-y[1:50000,]
  #   }
  #   y1<-y
  #   if (input$transpose1) {
  #     #message("Transposing Exp1...")
  #     yn=y[,1]
  #     cn=colnames(y1)[-1]
  #     y<-as.data.frame(t(y[,2:ncol(y)]))
  #     y<-setNames(y,yn)
  #     y<-cbind(cn,y)
  #     colnames(y)[1]<-"Pat"
  #   }
  #   #if (input$search > oldsearch)
  #   #{
  #   #	#message("New Search!")
  #   #   }
  #   #  oldsearch<<-input$search
  #   return(y)})
  # evidInput<-reactive({
  #   #message("READING EVIDENCE")
  #   evidFile <- input$evidence
  #   
  #   if (is.null(evidFile))
  #     return(NULL)
  #   z<-read.csv(mrnaFile$datapath, check.names=FALSE,header=input$header, sep=input$sep, quote=input$quote)
  #   return(z)})
  # mrnaInput<-reactive({
  #   #message("READING MRNAINPUT")
  #   mrnaFile <- input$mrnafile
  #   
  #   if (!is.null(mrnaFile))
  #     z<-read.csv(mrnaFile$datapath, check.names=FALSE,header=input$header, sep=input$sep, quote=input$quote)
  #   else if (input$loaddemo > 0) {
  #     #message("DEFAULT DATASET")
  #     z<-read.csv("./data/peax_demo_mir.csv", header=TRUE, sep=",", quote="\"")
  #   }
  #   else
  #     return(NULL)
  #   if (input$transpose2) {
  #     #message("Transposing Exp2...")
  #     zn=z[,1]
  #     z<-as.data.frame(t(z[,2:ncol(z)]))
  #     z<-setNames(z,zn)
  #   }
  #   return(z)})
  
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  ################
  get_splitpheno<-function(x,node,tr)
  {
    #message(paste0("DEBUG:get_splitpheno",node," tree:",tr))
    retpheno<-x
    #
    #  1
    #2  3
    #45  67
    curvar<-paste0("var",tr,"_",node%/%2)
    cursld<-paste0("range_slider",tr,"_",node%/%2)
    
    if (node==1) {
      return(x)
    }
    if (is.null(input[[curvar]]) || is.null(input[[cursld]])) {
      #message("NULL HERE")
      return(x)
    }
    slval<-input[[cursld]]
    #x
    
    if (node>1) {
      retpheno<-get_splitpheno(x,node%/%2,tr) #get parent pheno
      if (is.null(retpheno) || nrow(retpheno)==0) {
        retpheno<-retpheno
      }
      else if (node%%2==0) { #odd/left branch
        #message("Debug left")
        retpheno<-retpheno[retpheno[input[[curvar]]]<slval,]
        #message("Debug leftdone")
      }
      else {
        #message("Debug right")
        retpheno<-retpheno[retpheno[input[[curvar]]]>=slval,]
        #message("Debug rtdone")
      }
    }
    else
    {
      retpheno<-x  ## everybody
    }
    ds<-input[[paste0("depth_slider",tr)]]
    #  rr$Group=4
    # if it's a leaf, add the group annotation
    leaf_thresh=2^ds
    
    if (!is.null(retpheno) && nrow(retpheno)>0) {
      retpheno$Group=node%%leaf_thresh+1
    }
    else {
      #message("DEBUG:returning NULL retpheno")
    }
    return(retpheno)
  }
  
  # getPhenoExp<-function(x, curtr) {
  #   
  #   y <- mirInput()
  #   
  #   ds1<-input$depth_slider1
  #   #ds2<-input$depth_slider2
  #   
  #   #if (curtr==1) {
  #     oldselmir=oldselmir1
  #     ds<-input$depth_slider1
  #     tbl<-input$testTbl1_1
  #   #}
  #   # else
  #   # {
  #   #   oldselmir=oldselmir2
  #   #   ds<-input$depth_slider2
  #   #   tbl<-input$testTbl2_1
  #   # }
  #   
  #   #selmir<-tbl[curtr]
  #   selmir<-tbl
  #   # BROKEN here? selmir<-"rand.miR.101"
  #   #message("getPhenoExp")
  #   #message(selmir)
  #   #message("getPhenoExp2")
  #   if (is.null(selmir) || is.na(selmir)) {
  #     if (is.null(oldselmir) || oldselmir==0)
  #       selmir = colnames(y)[2]
  #     else
  #       selmir = oldselmir
  #   }
  #   else {
  #     #selmir=colnames(y)[selmir]
  #     selmir=selmir
  #   }
  #   
  #   #make a simple pheno vector with PT column and pheno value
  #   themir<-cbind(y[1],y[selmir])
  #   colnames(themir)[2]<-"Pheno"
  #   
  #   #message("DEBUG:getPheno3d")
  #   x<-merge(x,themir,all.x=TRUE,by=1) #merge pheno by Pt # in first col
  #   #message("DEBUG:getPhenoExp3d")
  #   return(x)
  # }
  
  
  getPheno<-function(curtr,histvar) {
    #message(paste0("DEBUG:getPheno tree:",curtr))
    #ds1<-input$depth_slider1
    #ds2<-input$depth_slider2
    
    #if (curtr==1) {
      oldselmir=oldselmir1
      ds<-input$depth_slider1
      tbl<-input$testTbl1_1
    # }
    # else
    # {
    #   oldselmir=oldselmir2
    #   ds<-input$depth_slider2
    #   tbl<-input$testTbl2_1
    # }
    
    if (is.null(ds)) {
      return(NULL)
    }
    
    x <- datasetInput()
    #message("DEBUG:getPheno2")
    if (is.null(x)) {
      #message("DEBUG:datasetInput is null")
      return(NULL)
    }
    
    #x<-x[c(colnames(x)[1],sapply(1:(2^ds-1),function(x) input[[paste0("var",curtr,"_",x)]]))]
    x<-x[c(input[[paste0("hist_input",histvar)]],colnames(x)[1],sapply(1:(2^ds-1),function(x) input[[paste0("var",curtr,"_",x)]]))]
    #x<-x[c(input[["hist_input2"]],colnames(x)[1],sapply(1:(2^ds-1),function(x) input[[paste0("var",curtr,"_",x)]]))]
    x <- x[complete.cases(x), ]
    #message("DEBUG:getPheno3")
    #y <- mirInput()
    
    # bind all leaves
    #e.g. depth 3
    # 8 leaf nodes, labeled 8 through 15
    
    leaf_nodes<-2^(ds)
    ret<-NULL
    #message("DEBUG:getPheno3e")
    for (lf in leaf_nodes:(leaf_nodes*2-1)) {
      ret<-rbind(ret,get_splitpheno(x,lf,curtr))
    }
    #message("DEBUG:getPheno3f")
    if (!is.null(ret)) {
      ret<-ret[with(ret,order(ret[,1])),]
    }
    #message("DEBUG:getPheno3g")
    #message("DEBUG:Leaving getPheno")
    #message(paste0("GetPheno: ",ret))
    return(ret)
  }
  
  # getPheno1<-reactive({
  #   return(getPheno(1))
  # })
  # 
  # getPheno2<-reactive({
  #   return(getPheno(2))
  # })
  
  
  
  # testcor <-function(curtr)
  # {
  #   #if (tr==1) {
  #   #  data<-getPheno1()
  #   #}
  #   #else {
  #   #  data<-getPheno2()
  #   #}
  #   #if (is.null(data)) {
  #   #	return(NULL)
  #   #}
  #   #message("DEBUG:Enter testcor()...")
  #   z<-mrnaInput()
  #   y<-mirInput()
  #   
  #   if (is.null(z) || is.null(y)) {
  #     return(NULL)
  #   }
  #   y<-y[y[,1] %in% z[,1],]
  #   z<-z[z[,1] %in% y[,1],]
  #   
  #   # do we need Group for corr?
  #   #z$Group<-data$Group
  #   
  #   # remove patients (col 1) before calculating
  #   
  #   if (curtr==1)
  #     tbl<-input$testTbl1_1
  #   else
  #     tbl<-input$testTbl2_1
  #   
  #   selmir<-tbl[curtr]
  #   
  #   if (is.null(selmir) || is.na(selmir))
  #   {
  #     #message("selmir undefined in testcor")
  #     selmir<-oldselmir
  #   }
  #   #message(paste0("DEBUG:Starting correlation...:",selmir))
  #   pvals<-cor(y[selmir],z[-1],method="spearman")
  #   #message("DEBUG:Ending correlation...")
  #   #pvals<-pvals[,abs(pvals)>0.6]
  #   pvals<-sort(pvals[1,])
  #   #nm<-names(pvals)
  #   nm<-names(pvals)
  #   zcor=data.frame(Probe = nm,Spearman=round(as.numeric(pvals),4))
  #   #message("DEBUG:Leaving testcor()...")
  #   return(zcor)
  # }
  
  
  # testdf1 <- reactive(
  #   {
  #     return(testdf(1))
  #   })
  # testdf2 <- reactive(
  #   {
  #     return(testdf(2))
  #   })
  # 
  # testdf<-function(tr) {
  #   ts<-input$tabSelected
  #   #message(ts)
  #   a <- input$range_slider1_1
  #   
  #   if (tr==1) {
  #     data<-getPheno1()
  #   }
  #   else {
  #     data<-getPheno2()
  #   }
  #   if (is.null(data)) {
  #     #message("DEBUG:NULL data in getPheno")
  #     return(NULL)
  #   }
  #   #message(paste0("DEBUG:Enter testdf()... with tree:",tr))
  #   y<-mirInput()
  #   data<-data[data[,1] %in% y[,1],]
  #   y<-y[y[,1] %in% data[,1],]
  #   
  #   y$Group<-data$Group
  #   y$Group<-as.factor(y$Group)
  #   yaov<-0
  #   
  #   #message("DEBUG:begin Anova step")
  #   
  #   y$Group<-data$Group
  #   y<-data.frame(as.matrix(y[-1]))
  #   y$Group<-as.factor(y$Group)
  #   Klist<-sort(unique(y$Group))
  #   #message(paste0("DEBUG:Klist=",Klist))
  #   K<-length(Klist)
  #   ### FILTER HERE
  #   Yib<-mclapply(Klist,function(tmp) colMedians(y[y$Group==tmp,][-ncol(y)]))
  #   
  #   if (input$searchall > oldsearchall) {
  #     oldsearchall<<-input$searchall
  #     oldsearch<<-input$search
  #     curfilt<<-NULL
  #   }
  #   
  #   
  #   if (input$search > oldsearch)
  #   {
  #     #message(paste0("Search:",input$search," ",oldsearch))
  #     oldsearch<<-input$search
  #     curfilt<<-input$testTbl1_1
  #     if (is.null(curfilt)) {
  #       curfilt<<-oldselmir1
  #     }
  #     #message(paste0("New Search!",curfilt))
  #   }
  #   if (!is.null(curfilt)) {
  #     Ydf<-data.frame(Yib)
  #     curthresh<-input$fold_thresh/100
  #     Yfilt<-mclapply(Klist,function(tmp) median(y[y$Group==tmp,curfilt]))
  #     Yfilt<-(Yfilt<(1-curthresh))*-1+(Yfilt>(1+curthresh)*1)
  #     Ydf<-data.frame((Ydf<(1-curthresh))*-1+(Ydf>(1+curthresh)*1))
  #     found<-t(data.frame(colSums(t(Ydf)==Yfilt)))
  #     names(found)<-colnames(y)[-ncol(y)]
  #     found<-found[found==K]
  #     found<-names(found)
  #     y<-y[found]
  #     y$Group<-as.factor(data$Group)
  #     
  #   }
  #   else {
  #   }
  #   
  #   nm<-names(y)[-ncol(y)]
  #   #nm<-nm[-1]
  #   N<-nrow(y)
  #   Ybar<-colMeans(y[-ncol(y)])
  #   Yib<-mclapply(Klist,function(tmp) colMeans(y[y$Group==tmp,][-ncol(y)]))
  #   ni<-mclapply(Klist,function(tmp) nrow(y[y$Group==tmp,]))
  #   expV<-mclapply(1:K,function(tmp) (ni[[tmp]]*(Yib[[tmp]]-Ybar)^2)/(K-1))
  #   expVar<-0
  #   for (i in 1:length(expV)) {
  #     expVar=expVar+expV[[i]]
  #   }
  #   
  #   #yy<-mclapply(Klist,function(tmp) apply(y[y$Group==tmp,][-ncol(y)],2,'-',Yib[[tmp]]))
  #   unexpV<-0
  #   for (i in 1:K) {
  #     yy<-y[y$Group==Klist[i],][-ncol(y)]-rep.row(Yib[[i]],ni[[i]])
  #     yy<-yy^2/(N-K)
  #     unexpV=unexpV+colSums(yy)
  #   }
  #   
  #   Ft<-expVar/unexpV
  #   
  #   pval<-(1-pf(Ft,df1=K-1,df2=N-K))
  #   #message("DEBUG:end Anova step")
  #   numtests<<-numtests+length(y)
  #   
  #   pvals<-pval
  #   yaov=data.frame(Probe = nm,p_val=round(as.numeric(pvals),5))
  #   #WHY REMOVE THE LAST ROW?yaov=yaov[-1,]
  #   yaov=yaov[with(yaov,order(p_val)),]
  #   #message("DEBUG:Leaving testdf()...")
  #   return(yaov)
  # }
  
  # numTests <- reactive({
  #   a <- input$range_slider1_1
  #   b <- input$range_slider1_2
  #   c <- input$range_slider1_3
  #   d <- input$range_slider1_1
  #   e <- input$range_slider1_2
  #   f <- input$range_slider1_3
  #   g <- input$searchall
  #   h <- input$search
  #   paste0("#Tests:",numtests)
  # })
  # 
  # output$numtests<-renderText({
  #   numTests()
  # })
  
  for (tr in 1:1) {
    local ({
      curtr<-tr
      # tbl1<-paste0("testTbl",curtr,"_1")
      # #tbl2<-paste0("testTbl",curtr,"_2")
      # # messages MIR table with selectable rows and conditional formatting
      # output[[tbl1]] <- renderUI({
      #   tdf<-testdf(curtr)
      #   if (!is.null(tdf)) {
      #     # restrict to top 32 mirs
      #     tdf<-na.omit(tdf[1:32,])
      #     tdf[,3]<-tdf[,2]
      #     colnames(tdf)[3]<-colnames(tdf)[2]
      #     colnames(tdf)[2]<-"Link"
      #     mbnames<-tdf
      #     mbnames[,1]<-gsub("_st","",mbnames[,1])
      #     mbnames[,1]<-gsub("[.]","_",mbnames[,1])
      #     mbnames[,1]<-gsub("_star","_5p",mbnames[,1])
      #     tdf[,2]<-paste("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",mbnames[,1]," target=_blank>",mbnames[,1],"</a>",sep="")
      #     HTML(df2html(tdf, class = "tbl selRow", id = tbl1,
      #                  cellClass = cbind(rep(NA, nrow(tdf)), rep(NA, nrow(tdf)), ifelse(tdf[,3]>=PTHRESH[input$pv_thresh], 'cellRed', 'cellGreen'))
      #     )
      #     )
      #   }
      # })
      # HERE WE message OUT MICRORNA
      # messages table with selectable rows and conditional formatting
      # output[[tbl2]] <- renderUI({
      #   tc<-NULL
      #   
      #   if (curtr==1)
      #     tbl<-input$testTbl1_1
      #   else
      #     tbl<-input$testTbl2_1
      #   selmir<-tbl[curtr]
      #   if (!is.null(tbl1))
      #     tc<-testcor(curtr)
      #   if (!is.null(tc)) {
      #     # restrict to top 10 mRNA
      #     tc<-na.omit(tc[1:10,])
      #     tc[,3]<-tc[,2]
      #     #tc[,4]<-tc[,2]
      #     colnames(tc)[3]<-colnames(tc)[2]
      #     colnames(tc)[2]<-"Link"
      #     
      #     mbnames<-tc
      #     mbnames[,1]<-gsub("_st","",mbnames[,1])
      #     mbnames[,1]<-gsub("[.]","_",mbnames[,1])
      #     mbnames[,1]<-gsub("_star","_5p",mbnames[,1])
      #     tc[,2]<-paste("<a href=http://www.mirbase.org/cgi-bin/query.pl?terms=",mbnames[,1]," target=_blank>","MiRBase","</a>",sep="")
      #     mbnames[,1]<-gsub("mir","miR",mbnames[,1])
      #     mbnames[,1]<-gsub("_","-",mbnames[,1])
      #     #tc[,4]<-sapply(1:nrow(tc), function(x) {
      #     #cur_mm<-get.multimir(target=oldselmir1,mirna=mbnames[x,1])$validated$pubmedid;
      #     #	cur_mm<-""
      #     #paste("<a href=www.ncbi.nlm.nih.gov/pubmed/",cur_mm[1],">",length(cur_mm),"</a>",sep="")
      #     
      #     #})
      #     #colnames(tc)[4]<-"Bind"
      #     HTML(df2html(tc, class = "tbl selRow", id = tbl2,
      #                  cellClass = cbind(rep(NA, nrow(tc)), rep(NA, nrow(tc)), ifelse(abs(tc[,3])>=input$sp_thresh, 'cellGreen', 'cellRed'),NA)
      #                  #          cellClass = cbind(rep(NA, nrow(tc)), rep(NA, nrow(tc)), ifelse(abs(tc[,3])>=input$sp_thresh, 'cellGreen', 'cellRed'))
      #     )
      #     )
      #   }
      # })
      # 
      # curchoose_mir<-paste0("choose_mir",curtr)
      # output[[curchoose_mir]] <- renderUI({
      #   if (is.null(input$mirfile))
      #     return(NULL)
      #   
      #   selectInput("mir", "microRNA", colnames(mirInput()))
      # })
      # 
      # # Compute the forumla text in a reactive expression since it is
      # formulaText1 <- reactive({
      #   selmir<-input$testTbl1_1
      #   if (selmir==null) {
      #     selmir = y[1,1]
      #   }
      #   paste(selmir," ~","pheno")
      # })
      # # Compute the forumla text in a reactive expression since it is
      # formulaText2 <- reactive({
      #   selmir<-input$testTbl2_1
      #   if (selmir==null) {
      #     selmir = y[1,1]
      #   }
      #   paste(selmir," ~","pheno")
      # })
      
      ds1<-2
      # Generate a plot of the requested variable against x and only
      # include outliers if requested
      for (i in 1:2^ds1) {
        local ({
          local_i <- i
          curPlot<-paste0("plotVar1",local_i)
          output[[curPlot]] <- renderPlot({
            message("DEBUG1: render hist Plot")
            selVariable <- input[["hist_input1"]]
            #message(paste0("DEBUG1: selVariable  ",selVariable))
            #bp<-getPhenoExp(getPheno(1,1),1)
            bp<-getPheno(1,1)
            #message(paste0("DEBUG1: getPheno(1)  ",getPheno(1)))
            #message(paste0("DEBUG1: bp  ",bp))
            bp<-bp[bp$Group==local_i,]
            #write.csv(bp, file = paste0("Cohort",local_i),row.names = FALSE)
            bp<-bp[selVariable]
            if (length(bp)>0) {
              hist(bp[,1],main=paste0("Distribution of ",selVariable),xlab=selVariable,axes=TRUE,right=FALSE)
            }
          },width=200,height=250) # this is the actual size of the plot output
        })
      }
      
      # Generate a plot of the requested variable against x and only
      # include outliers if requested
      for (i in 1:2^ds1) {
        local ({
          local_i <- i
          curPlot<-paste0("plotVar2",local_i)
          output[[curPlot]] <- renderPlot({
            message("DEBUG2: render hist Plot")
            selVariable2 <- input[["hist_input2"]]
            #message(paste0("DEBUG2: selVariable  ",selVariable2))
            #bp<-getPhenoExp(getPheno(1,2),1)
            bp<-getPheno(1,2)
            #message(paste0("DEBUG2: getPheno(1)  ",getPheno(1)))
            #message(paste0("DEBUG2: bp  ",bp))
            bp<-bp[bp$Group==local_i,]
            #message(paste0("DEBUG2: selVariable BEFORE ",selVariable2))
            #message(paste0("DEBUG2: selVariable BEFORE ",colnames(bp)))
            bp<-bp[selVariable2]
            #message(paste0("DEBUG2: selVariable AFTER ",selVariable2))
            if (length(bp)>0) {
              hist(bp[,1],main=paste0("Distribution of ",selVariable2),xlab=selVariable2,axes=TRUE,right=FALSE)
            }
          },width=200,height=250) # this is the actual size of the plot output
        })
      }
      
      #################
      
    })
  } # done with both trees
  
  # observe ({
  #   #message(paste0("Table 5: ",  ifelse(is.null(input$testTbl1_2), "NULL", input$testTbl1_2 )))
  #   mrna_sel=input$testTbl1_2
  #   #message(input$testTbl1_2)
  #   message(numtests)
  # })
  # observe ({
  #   if (!is.null(input$testTbl1_1))
  #     oldselmir1<<-input$testTbl1_1
  #   if (!is.null(input$testTbl2_1))
  #     oldselmir2<<-input$testTbl2_1
  #   
  # })
  
  # input$file1 will be NULL initially. After the user selects and uploads a
  # file, it will be a data frame with 'name', 'size', 'type', and 'datapath'
  # columns. The 'datapath' column will contain the local filenames where the
  # data can be found.
  
  output$pcontents <- renderTable({
    datasetInput()[c(1:80),]
  })
  # output$mircontents <- renderTable({
  #   mirInput()
  # })
  # output$mrnacontents <- renderTable({
  #   mrnaInput()
  # })
  
  
  
  
  messagesplit<-function(node,samples,varname,varval) {
    eqop="<"
    hdr=""  # should only be first split
    if (node%%2==1) { #right
      hdr=""
    }
    tmpstr<-paste0(hdr,"{\"samples\":",samples,",\"value\":[0],\"label\":\"","[",node,"]  ",varname,eqop,varval,"\",\"type\":\"split\",\"children\":[")
    return (tmpstr)
  }
  
  messageleaf<-function(node,samples,tr) {
    hdr=""
    clsr=""
    ds<-input[[paste0("depth_slider",tr)]]
    if ((node%%2)==0) { #right, because pheno labeling is 1-based
      hdr=""
      if (node>1) { #2,4,8
        clsrcnt<-0
        for (i in 1:ds) {
          clsrcnt<-clsrcnt+(node%%(2^i)==0) # add some for 2,4,8,16
        }
        for (i in 1:clsrcnt) {
          clsr<-paste0(clsr,"]}")
        }
      }
      else
        clsr<-paste0(clsr,"]}")
    }
    tmpstr<-paste0(hdr,"{\"samples\":",samples,",\"value\":[0],\"label\":\"","Pheno",node,"\"",",\"type\":\"leaf\"}",clsr)
    return (tmpstr)
  }
  
  #   1
  # 2   3
  #45   67
  #89abcdef
  
  #124895ab36cd7ef
  message_bst<-function(idx,sz) {
    tidx<-c(idx)
    if (idx*2<sz) {
      tidx<-c(tidx,message_bst(idx*2,sz))
      tidx<-c(tidx,message_bst(idx*2+1,sz))
    }
    return(tidx)
  }
  
  tree_traverse<-function(x,tr) {
    
    retstr<-""
    #x <- datasetInput()
    #message("HERE IN tree_traverse")
    ds<-input[[paste0("depth_slider",tr)]]
    if (is.null(ds)) {
      return(retstr)
    }
    x<-x[c(colnames(x)[1],sapply(1:(2^ds-1),function(nx) input[[paste0("var",tr,"_",nx)]]))]
    if (is.null(x)) {
      return(retstr)
    }
    x <- x[complete.cases(x), ]
    curpheno<-x
    
    numnodes=2^(ds+1)
    tidx<-message_bst(1,numnodes)
    leaf_thresh=2^(ds)
    
    i<-1
    while (i <(numnodes-1)) {
      idx<-tidx[i]
      curpheno<-get_splitpheno(x,idx,tr)
      
      if (idx>=leaf_thresh) {
        phenonum<-(idx%%leaf_thresh)+1
        # left leaf
        retstr<-paste0(retstr,messageleaf(phenonum,nrow(curpheno),tr),",")
        i=i+1
        idx<-tidx[i]
        phenonum<-(idx%%leaf_thresh)+1
        curpheno<-get_splitpheno(x,idx,tr)
        retstr<-paste0(retstr,messageleaf(phenonum,nrow(curpheno),tr))
        # right leaf
        if (i<(numnodes-1)) {
          # up and over
          retstr<-paste0(retstr,",")
          i=i+1
          idx<-tidx[i]
        }
      }
      else {
        retstr<-paste(retstr,messagesplit(idx,nrow(curpheno),input[[paste0("var",tr,"_",idx)]],input[[paste0("range_slider",tr,"_",idx)]]),sep="")
        # traverse down
        i=i+1
        idx<-tidx[i]
      }
    }
    
    message (retstr)
    return(retstr[1])
  }
  ########### end tree_traverse
  
  for (tr in 1:1) { # 2 trees
    for (i in 1:n) {
      #message(paste0("DEBUG: In outer loop",tr," ",i))
      local({
        curtr<-tr
        #make dynamic slider
        row <- i
        #message(paste0("DEBUG: In loop",curtr," ",row))
        curvar<-paste0("var",curtr,'_',row)
        #message("MAH in loopB4")
        curcol<-paste0('choose_columns',curtr,'_',row)
        hist_var1<-'hist_var1'
        hist_input1<-'hist_input1'
        hist_var2<-'hist_var2'
        hist_input2<-'hist_input2'
        output[[curcol]] <- renderUI({
          #message("DEBUG: render cur col")
          if (is.null(input$file1) && (input$loaddemo==0))
          {#message("MAH in null nputfile")
            return(NULL)
          }
          selectInput(curvar, paste0("Variable ",row), colnames(datasetInput()))
        })
        output[[hist_var1]] <- renderUI({
          #message("DEBUG: render cur col")
          if (is.null(input$file1) && (input$loaddemo==0))
          {#message("MAH in null nputfile")
            return(NULL)
          }
          selectInput(hist_input1, ("Histogram Var 1"), colnames(datasetInput()))
        })
        output[[hist_var2]] <- renderUI({
          #message("DEBUG: render cur col")
          if (is.null(input$file1) && (input$loaddemo==0))
          {#message("MAH in null nputfile")
            return(NULL)
          }
          selectInput(hist_input2, ("Histogram Var 2"), colnames(datasetInput()))
        })
        id<-paste0('range_slider',curtr,'_',row)
        output[[id]] <- renderUI({
          #message("DEBUG: render depth slider")
          depth<-input[[paste0("depth_slider",curtr)]]
          if (is.null(input$file1) && (input$loaddemo==0))
            return(NULL)
          if (is.null(input[[curvar]]))
            return(NULL)
          x <- datasetInput()
          x<-x[c(colnames(x)[1],sapply(1:(2^depth-1),function(nx) input[[paste0("var",curtr,"_",nx)]]))]
          #message(paste0("DEBUG:size of x is",nrow(x)))
          x <- x[complete.cases(x), ]
          #message(paste0("DEBUG:complete cases size of x is",nrow(x)))
          splitx<-get_splitpheno(x,row,curtr)
          if (!is.null(splitx) && (nrow(splitx)>0)) {
            splitx<-splitx[input[[curvar]]]
            slmin1 <- floor(min(splitx[input[[curvar]]]))
            slmax1 <- ceiling(max(splitx[input[[curvar]]]))
            name_idx=grep(input[[curvar]],colnames(splitx),fixed=TRUE)
            slmed1 <- ceiling(as.numeric(sapply(splitx[name_idx],median)))
          }
          else {
            #message(paste0("DEBUG: Null when creating slider:",id))
            slmin1<-0
            slmax1<-0
            slmed1<-0
          }
          sliderInput(inputId = id,
                      label = paste(""),
                      #min = slmin1+(row<=(2^(depth-2)+1)), max = slmax1+(row>(2^(depth-2)+1)), value = slmed1)
                      min = slmin1, max = slmax1+1, value = slmed1)
        })
        curhist<-paste0("hist",curtr,"_",row)
        
        output[[curhist]] <- renderPlot({
          #message("DEBUG: render hist Plot")
          x <- datasetInput()
          depth<-input[[paste0("depth_slider",curtr)]]
          x<-x[c(colnames(x)[1],sapply(1:(2^depth-1),function(nx) input[[paste0("var",curtr,"_",nx)]]))]
          x <- x[complete.cases(x), ]
          splitx<-get_splitpheno(x,row,curtr)
          curnode<-splitx
          if (is.null(curnode) || nrow(curnode)==0) {
            #message("DEBUG: Null histogram")
          } else {
            bp<-curnode[input[[curvar]]]
            if (length(bp)>0) {
              par(pin=c(0.3,0.3))
              par(mar=c(0,0,0,0))
              hist(bp[,1],main=NULL,xlab=NULL,ylab=NULL,axes=FALSE,right=FALSE)
            }
          }
        },width=200,height=20)
      })
    }
  }
  
  
  
  observe({
    #message("DEBUG:range_slider observe")
    
    for (t in 1:1) {
      for (i in 1:n) {
        input[[paste0("range_slider",t,"_",i)]]
      }}
    for (i in 1:n) {
      #message(input[[paste0("range_slider1_",i)]])
    }
    input$range_slider1_1 # Do take a dependency on input$saveButton
    input$range_slider1_2 # Do take a dependency on input$saveButton
    input$range_slider1_3 # Do take a dependency on input$saveButton
    
    
    # isolate a whole block
    data <- isolate({
      #data <- ({
      lvls<-input$depth_slider1
      lvls2<-input$depth_slider2
      #message("IN data block")
      a <- input$range_slider1_1
      b <- input$range_slider1_2
      c <- input$range_slider1_3
      d <- input$range_slider2_1
      e <- input$range_slider2_2
      f <- input$range_slider2_3
      x <- datasetInput()
      
      
      pattern<-tree_traverse(x,1)
      sink("./www/tmp.json")
      cat(pattern)
      sink()
      
      file.rename("./www/tmp.json",paste0("./www/tree",thetime,"_1.json"))
      
      pattern<-tree_traverse(x,2)
      sink("./www/tmp.json")
      cat(pattern)
      sink()
      file.rename("./www/tmp.json",paste0("./www/tree",thetime,"_2.json"))
      
    })
  })
})
