run.gui <-
function(){

    e1 <- new.env(parent = .GlobalEnv)  
    e2 <- new.env(parent = .GlobalEnv)  
    e3 <- new.env(parent = .GlobalEnv)  

    base <- tktoplevel()
    tkwm.title(base,"iDEMO")
    Mainfrm <- tkframe(base,relief="groove")
    mainfrm <- tkframe(Mainfrm,relief="groove",borderwidth=3)

    F <- tkframe(mainfrm,relief="groove",borderwidth=3)
    f0 <- tkframe(F)
    tkgrid(tklabel(f0,text="Linear Degradation Model",font=tkfont.create(family="times",size=18)), sticky="w")

    tkgrid(f0)

    nb <- tk2notebook(F, tabs = c("Basic information", "Parameter estimation", "Lifetime information"))

    tb1 <- tk2notetab(nb, "Basic information")

    f1 <- tkframe(tb1)
    datapick <- tclVar('<no data>')

    pickone <- function(){
        namelist <- ls(envir=e1)
        tt<-tktoplevel()
        tkwm.title(tt,"Name list")
        scr <- tkscrollbar(tt, repeatinterval=5,
                   command=function(...)tkyview(tl,...))
        tl<-tklistbox(tt,width=27,height=10,selectmode="single",yscrollcommand=function(...)tkset(scr,...),background="white")
        if(length(namelist)  == 0) {
            tkinsert(tl,"end", c())
        } else {
                for (i in (1:length(namelist))) {
                    tkinsert(tl,"end",namelist[i])
                }
        }

        tkselection.set(tl,0)

        OnOK <- function()
        {
		if(length(namelist)==0){
	            tclvalue(datapick) <- "<no data>"
		}else{
	            tclvalue(datapick) <- namelist[as.numeric(tkcurselection(tl))+1]
		}
            tkdestroy(tt)
        }
        OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
        tkgrid(tklabel(tt,text="Data sets"))
        tkgrid(tl,scr)
        tkgrid.configure(scr,rowspan=4,sticky="nsw")
        tkgrid(OK.but)
        tkfocus(tt)
    }

    data.lab <- tkbutton(f1,textvariable=datapick,relief='groove',fg="blue",
                  command=pickone)

    open <- function(){
        name <- tclvalue(tkgetOpenFile(
            filetypes = "{{Text Files} {.txt}} {{All files} *}"))
        if (name == "") return()
        Data <- read.table(name, header = TRUE)
	  if( any( is.na(Data) ) ) stop("Please check your data. Data with missing value is not available!")

        totalnumfile <- unlist(strsplit(name, "/"))
        endnumfile <- length(totalnumfile)
        filenamechar <- totalnumfile[endnumfile]

        if(length(Data)!=0){
	        tclvalue(datapick) <- substr(filenamechar, 1, nchar(filenamechar)-4)
	        if (!is.null(Data)) assign(tclvalue(datapick), Data, envir=e1)
	        if (!is.null(Data)) assign(tclvalue(datapick), Data, envir=.GlobalEnv)
	 }
    }
    Open.but <- tkbutton(f1,text=" Import data ",
      command=open
    )

    lab1 <- tklabel(f1,text="")
    tkgrid(lab1, sticky="w")


    f1.1 <- tkframe(tb1)
    plot.path <- function(){
		if(!exists(tclvalue(datapick))){
			stop("No data set! Please choose one to analyze!")
		}
		plot.degrad.path(tclvalue(datapick),get(tclvalue(datapick)))
    }
    plot.path.but <- tkbutton(f1.1,text="Show the degradation paths ", command=plot.path)

    box <- function(){
		if(!exists(tclvalue(datapick))){
			stop("No data set! Please choose one to analyze!")
		}
		idemo.box(tclvalue(datapick),get(tclvalue(datapick)))
    }
    box.but <- tkbutton(f1.1,text="Degradation data summary", command=box)

    pseudoest <- function(){
		if(!exists(tclvalue(datapick))){
			stop("No data set! Please choose one to analyze!")
		}
		if(!is.numeric(eval(parse(text=tclvalue(w)))) || length(eval(parse(text=tclvalue(w))))>1 || eval(parse(text=tclvalue(w)))<=0 ) stop("The inputted threshold must be numeric and a positive value, not a vector!")

		PseudoFailureTime(get(tclvalue(datapick)),eval(parse(text=tclvalue(w))),tclvalue(datapick))
    }
    pseudoest.but <- tkbutton(f1.1,text="   Pseudo failure time estimation   ", command=pseudoest)

    lab4 <- tklabel(f1.1,text='Threshold (\u03c9)')
    w <- tclVar("10 ")
    entry.Name4 <-tkentry(f1.1,width="25",textvariable=w)

    int.set <- function(){
	  if(!exists(tclvalue(datapick))){
		  stop("No data set! Please choose one to analyze!")
	  }
	    tt<-tktoplevel()
	    tkwm.title(tt,"Initial setting")
	    tt1 <- tkframe(tt,relief="groove",borderwidth=3)
	    cb0 <- tklabel(tt1,text='Drift rate (\u03b7)')
	    eta.int.wiz <- tclVar("0 ")
	    entry.Name0 <-tkentry(tt1,width="25",textvariable=eta.int.wiz)

	    lab00 <- tklabel(tt1,text="Variation sources",fg="brown", state='normal')

	    cb1 <- tklabel(tt1,text='Unit-to-unit variation (\u03c3_\u03b7)')
	    seta.int.wiz <- tclVar("1 ")
	    entry.Name1 <-tkentry(tt1,width="25",textvariable=seta.int.wiz)

	    cb2 <- tklabel(tt1,text='Time-dependent structure (\u03c3_B)    ')
	    sb.int.wiz <- tclVar("1 ")
	    entry.Name2 <-tkentry(tt1,width="25",textvariable=sb.int.wiz)

	    cb3 <- tklabel(tt1,text='Measurement error (\u03c3_\u03b5)')
	    se.int.wiz <- tclVar("1 ")
	    entry.Name3 <-tkentry(tt1,width="25",textvariable=se.int.wiz)

	    rb1 <- tkradiobutton(tt1)
	    rb2 <- tkradiobutton(tt1)
	    rb3 <- tkradiobutton(tt1)
	    rb4 <- tkradiobutton(tt1)
	    rb5 <- tkradiobutton(tt1)
	    rbValue <- tclVar("Nelder-Mead")
	    tkconfigure(rb1,variable=rbValue,value="Nelder-Mead",text="Nelder-Mead")
	    tkconfigure(rb2,variable=rbValue,value="BFGS",text="BFGS")
	    tkconfigure(rb3,variable=rbValue,value="CG",text="CG")
	    tkconfigure(rb4,variable=rbValue,value="L-BFGS-B",text="L-BFGS-B")
	    tkconfigure(rb5,variable=rbValue,value="SANN",text="SANN")

	    tkgrid(lab01 <- tklabel(tt1,text="Model Parameters",fg="blue"),lab02 <- tklabel(tt1,text="Initial value"), sticky="ew")
	    tkgrid.configure(lab01,sticky="w")
  	    tkgrid(cb0, entry.Name0, sticky="w")
	    tkgrid(lab00, sticky="w")
	    tkgrid(cb1, entry.Name1, sticky="w")
	    tkgrid(cb2, entry.Name2, sticky="w")
	    tkgrid(cb3, entry.Name3, sticky="w")
	    tkgrid(lab30 <- tklabel(tt1,text="Optimization algorithm"), rb1, sticky="w")
	    tkgrid(tklabel(tt1,text=""), rb2, sticky="w")
	    tkgrid(tklabel(tt1,text=""), rb3, sticky="w")
	    tkgrid(tklabel(tt1,text=""), rb4, sticky="w")
	    tkgrid(tklabel(tt1,text=""), rb5, sticky="w")
	    tkgrid(tt1)


	    OnOK <- function()
	    {
		  assign('etaVal', eval(parse(text=tclvalue(eta.int.wiz))), envir=e3)
		  assign('setaVal', eval(parse(text=tclvalue(seta.int.wiz))), envir=e3)
		  assign('sbVal', eval(parse(text=tclvalue(sb.int.wiz))), envir=e3)
		  assign('seVal', eval(parse(text=tclvalue(se.int.wiz))), envir=e3)
		  assign('medthodChoice', as.character(tclvalue(rbValue)), envir=e3)

		  tkdestroy(tt)
		  Wiz()

	    }
	    OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
	    tkgrid(OK.but)
	    tkfocus(tt)

    }

    f1.2 <- tkframe(tb1)

    Wiz <- function() {
	  etaVal <- get('etaVal',envir=e3)
	  setaVal <- get('setaVal',envir=e3)
	  sbVal <- get('sbVal',envir=e3)
	  seVal <- get('seVal',envir=e3)
	  medthodChoice <- get('medthodChoice',envir=e3)

	  if(!is.numeric(etaVal) || length(etaVal)>1) stop("The inputted initial value must be numeric and a scale, not a vector!")
	  if(!is.numeric(setaVal) || length(setaVal)>1 || setaVal<0 ) stop("The inputted initial value must be numeric and a positive value, not a vector!")
	  if(!is.numeric(sbVal) || length(sbVal)>1 || sbVal<0 ) stop("The inputted initial value must be numeric and a positive value, not a vector!")
	  if(!is.numeric(seVal) || length(seVal)>1 || seVal<0 ) stop("The inputted initial value must be numeric and a positive value, not a vector!")

	  cat("Data:",tclvalue(datapick),"\n\n")
	  cat("Please wait...","\n\n")
        dist.wiz(get(tclvalue(datapick)),etaVal,setaVal,sbVal,seVal,medthodChoice)
	  if(length(medthodChoice)==0){
	        cat("Optimization Algorithm: Nelder-Mead","\n\n")
	  }else{
	        cat("Optimization Algorithm: ",medthodChoice,"\n\n")
	  }
    }
    wiz.cb <- tkcheckbutton(f1.2)
    sdma <- tklabel(f1.2,text="Single degradation model analysis",fg="blue",justify="left",state="disabled")
    sdma2 <- tklabel(f1.2,text="Specify a degradation model in the Parameter estimation tab  \nand choose  the  related product's  lifetime information from\nthe Lifetime information tab.",justify="left",state="disabled")
    cbValue00 <- tclVar("0")
    tkconfigure(wiz.cb,variable=cbValue00)
    wiz.but <- tkbutton(f1.1,text=" Degradation model selection ",command=int.set)

    tkgrid(tklabel(f1,text='Data set:'),data.lab,Open.but, sticky="w")
    tkgrid(tklabel(f1.1,text=''), sticky="w")
    tkgrid(tklabel(f1.1,text=''), plot.path.but, sticky="ew")
    tkgrid(tklabel(f1.1,text=''), sticky="w")
    tkgrid(tklabel(f1.1,text=''), box.but, sticky="ew")
    tkgrid(tklabel(f1.1,text=''), tklabel(f1.1,text=''), tklabel(f1.1,text='  '), lab4, sticky="ew")
    tkgrid(tklabel(f1.1,text=''), pseudoest.but, tklabel(f1.1,text='  '), entry.Name4, sticky="ew")
    tkgrid(tklabel(f1.1,text=''), sticky="w")
    tkgrid(tklabel(f1.1,text=''), wiz.but, sticky="ew")
    tkgrid(tklabel(f1.1,text=''), sticky="w")
    tkgrid(tklabel(f1.1,text=''), sticky="w")
    tkgrid(wiz.cb, sdma, sticky="w")
    tkgrid(tklabel(f1.2,text=''), sdma2, sticky="w")
    tkgrid.configure(wiz.cb,sticky="n")

    tb2 <- tk2notetab(nb, "Parameter estimation")
    f2 <- tkframe(tb2)

    tkgrid(lab01 <- tklabel(f2,text="Model Parameters",fg="blue", state='disabled'), lab02 <- tklabel(f2,text="Initial value", state='disabled'), sticky="ew")
    tkgrid.configure(lab01,sticky="w")

    cb0 <- tklabel(f2,text='Drift rate (\u03b7)', state='disabled')
    cbValue0 <- tclVar("1")
    tkconfigure(cb0)
    eta.int <- tclVar("0 ")
    entry.Name0 <-tkentry(f2,width="25",textvariable=eta.int, state='disabled')

    lab00 <- tklabel(f2,text="Variation sources",fg="brown", state='disabled')
    lab30 <- tklabel(f2,text="Optimization algorithm", state='disabled')
    cb1 <- tkcheckbutton(f2,text='Unit-to-unit variation (\u03c3_\u03b7)', state='disabled')
    cbValue1 <- tclVar("1")
    tkconfigure(cb1,variable=cbValue1)
    seta.int <- tclVar("1 ")
    entry.Name1 <-tkentry(f2,width="25",textvariable=seta.int, state='disabled')


    cb2 <- tkcheckbutton(f2,text='Time-dependent structure (\u03c3_B)    ', state='disabled')
    cbValue2 <- tclVar("1")
    tkconfigure(cb2,variable=cbValue2)
    sb.int <- tclVar("1 ")
    entry.Name2 <-tkentry(f2,width="25",textvariable=sb.int, state='disabled')


    cb3 <- tkcheckbutton(f2,text='Measurement error (\u03c3_\u03b5)', state='disabled')
    cbValue3 <- tclVar("1")
    tkconfigure(cb3,variable=cbValue3)
    se.int <- tclVar("1 ")
    entry.Name3 <-tkentry(f2,width="25",textvariable=se.int, state='disabled')

    rb1 <- tkradiobutton(f2)
    rb2 <- tkradiobutton(f2)
    rb3 <- tkradiobutton(f2)
    rb4 <- tkradiobutton(f2)
    rb5 <- tkradiobutton(f2)
    rbValue <- tclVar("Nelder-Mead")
    tkconfigure(rb1,variable=rbValue,value="Nelder-Mead",text="Nelder-Mead",state="disabled")
    tkconfigure(rb2,variable=rbValue,value="BFGS",text="BFGS",state="disabled")
    tkconfigure(rb3,variable=rbValue,value="CG",text="CG",state="disabled")
    tkconfigure(rb4,variable=rbValue,value="L-BFGS-B",text="L-BFGS-B",state="disabled")
    tkconfigure(rb5,variable=rbValue,value="SANN",text="SANN",state="disabled")

    tkgrid(cb0, entry.Name0, sticky="w")
    tkgrid(lab00, sticky="w")
    tkgrid(cb1, entry.Name1, sticky="w")
    tkgrid(cb2, entry.Name2, sticky="w")
    tkgrid(cb3, entry.Name3, sticky="w")
    tkgrid(lab30, rb1, sticky="w")
    tkgrid(tklabel(f2,text=""), rb2, sticky="w")
    tkgrid(tklabel(f2,text=""), rb3, sticky="w")
    tkgrid(tklabel(f2,text=""), rb4, sticky="w")
    tkgrid(tklabel(f2,text=""), rb5, sticky="w")

    tb3 <- tk2notetab(nb, "Lifetime information")
    f3 <- tkframe(tb3)
    f3.1 <- tkframe(tb3)
    f3.2 <- tkframe(tb3)

    tkgrid(aa2 <- tklabel(f3,text="Lifetime Distribution",fg="blue", state='disabled'), sticky="w")

    cb5 <- tkcheckbutton(f3,text='Mean-time-to-failure (MTTF)', state='disabled')
    cbValue5 <- tclVar("1")
    tkconfigure(cb5,variable=cbValue5)
    h <- tclVar("1e-10 ")

    cb6 <- tkcheckbutton(f3,text='qth quantile', state='disabled')
    cbValue6 <- tclVar("1")
    tkconfigure(cb6,variable=cbValue6)
    q <- tclVar("0.5")
    entry.Name6 <-tkentry(f3,width="25",textvariable=q, state='disabled')

    cb7 <- tkcheckbutton(f3,text="Show the plot of CDF estimation  ", state='disabled')
    cbValue7 <- tclVar("1")
    tkconfigure(cb7,variable=cbValue7)

    cb7.2 <- tkcheckbutton(f3,text="Show the plot of PDF estimation  ", state='disabled')
    cbValue7.2 <- tclVar("1")
    tkconfigure(cb7.2,variable=cbValue7.2)

    cdft1 <- tclVar("2000")
    cdft2 <- tclVar("12000")
    cdft3 <- tclVar("50")
    entry.Name7.1 <-tkentry(f3.2,width="25",textvariable=cdft1, state='disabled')
    entry.Name7.2 <-tkentry(f3.2,width="25",textvariable=cdft2, state='disabled')
    entry.Name7.3 <-tkentry(f3.2,width="25",textvariable=cdft3, state='disabled')

    lab6.91 <- tklabel(f3.1,text="Searching interval for the quantile", fg="brown", state='disabled')
    lab6.92 <- tklabel(f3.1,text="(domain of the plot)", fg="brown", state='disabled')
    lab7.from <- tklabel(f3.2,text="From", state='disabled')
    lab7.to <- tklabel(f3.2,text="To", state='disabled')
    lab7.incr <- tklabel(f3.2,text="Increment", state='disabled')

    lab7.2 <- tklabel(f2,text="Significant level (\u03b1)", state='disabled')

    lab8 <- tklabel(f2,text='Confidence Interval (CI)',fg="blue", state='disabled')
    alpha <- tclVar("0.05")
    entry.Name8 <-tkentry(f2,width="25",textvariable=alpha, state='disabled')

    lab8.1 <- tklabel(f2,text='Information matrix type',fg="brown", state='disabled')

    cb9 <- tkcheckbutton(f2,text="Fisher's information matrix (FIM)  ", state='disabled')
    cbValue9 <- tclVar("1")
    tkconfigure(cb9,variable=cbValue9)

    lab9.1 <- tklabel(f2,text='Observed information matrix (OIM)',fg="brown", state='disabled')

    cb10 <- tkcheckbutton(f2,text='Hessian matrix', state='disabled')
    cbValue10 <- tclVar("0")
    tkconfigure(cb10,variable=cbValue10)

    cb11 <- tkcheckbutton(f2,text='Score vector', state='disabled')
    cbValue11 <- tclVar("0")
    tkconfigure(cb11,variable=cbValue11)

    cb12 <- tkcheckbutton(f2,text='Robust matrix', state='disabled')
    cbValue12 <- tclVar("0")
    tkconfigure(cb12,variable=cbValue12)

    tkgrid(cb7.2, sticky="w")
    tkgrid(cb7, sticky="w")
    tkgrid(cb5, sticky="w")
    tkgrid(cb6, entry.Name6, sticky="w")
    tkgrid(lab6.91, sticky="w")
    tkgrid(lab6.92, sticky="w")
    tkgrid(lab7.from, lab7.01 <- tklabel(f3.2,text="     "), entry.Name7.1, sticky="w")
    tkgrid(lab7.to, lab7.02 <- tklabel(f3.2,text="     "), entry.Name7.2, sticky="w")
    tkgrid(lab7.incr, lab7.03 <- tklabel(f3.2,text="     "), entry.Name7.3, sticky="w")
    tkgrid(tklabel(f2,text=''), sticky="w")
    tkgrid(lab8, sticky="w")
    tkgrid(lab7.2, entry.Name8, sticky="w")
    tkgrid(lab8.1, sticky="w")
    tkgrid(cb9, sticky="w")
    tkgrid(lab9.1, sticky="w")
    tkgrid(cb10, sticky="w")
    tkgrid(cb11, sticky="w")
    tkgrid(cb12, sticky="w")

    lab13 <- tklabel(f3.2,text='Goodness of Fit',fg="blue", state='disabled')

    cb14 <- tkcheckbutton(f3.2,text='Plot pseudo failure time (PFT)', state='disabled')
    cbValue14 <- tclVar("1")
    tkconfigure(cb14,variable=cbValue14)

    lab14 <- tklabel(f3.2,text='Graphical representation',fg="brown", state='disabled')
    lab15 <- tklabel(f3.2,text='Hypothesis testing',fg="brown", state='disabled')

    cb15_1 <- tkcheckbutton(f3.2,text='Probability-probability plot', state='disabled')
    cbValue15_1 <- tclVar("0")
    tkconfigure(cb15_1,variable=cbValue15_1)

    cb15_2 <- tkcheckbutton(f3.2,text='Quantile-quantile plot', state='disabled')
    cbValue15_2 <- tclVar("0")
    tkconfigure(cb15_2,variable=cbValue15_2)

    cb15 <- tkcheckbutton(f3.2,text='Kolmogorov-Smirnov test', state='disabled')
    cbValue15 <- tclVar("1")
    tkconfigure(cb15,variable=cbValue15)


    cb16 <- tkcheckbutton(f3.2,text='Anderson-Darling test', state='disabled')
    cbValue16 <- tclVar("1")
    tkconfigure(cb16,variable=cbValue16)

    tkgrid(tklabel(f3.2,text=''), sticky="w")
    tkgrid(lab13, sticky="w")
    tkgrid(cb14, sticky="w")
    tkgrid(lab14, sticky="w")
    tkgrid(cb15_1, sticky="w")
    tkgrid(cb15_2, sticky="w")
    tkgrid(lab15, sticky="w")
    tkgrid(cb15, sticky="w")
    tkgrid(cb16, sticky="w")

    set.entry.state.wiz<-function(){
	    if(tclvalue(cbValue00)=="0"){
			tkconfigure(sdma,state = "normal")
			tkconfigure(sdma2,state = "normal")
			if(tclvalue(cbValue1)=="0" && tclvalue(cbValue2)=="0" && tclvalue(cbValue3)=="1"){
				tkconfigure(lab01,state = "normal")
				tkconfigure(lab02,state = "normal")
				tkconfigure(cb0,state = "normal")
				tkconfigure(entry.Name0,state = "normal")
				tkconfigure(lab00,state = "normal")
				tkconfigure(cb1,state = "normal")
				tkconfigure(entry.Name1,state = "disabled")
				tkconfigure(cb2,state = "normal")
				tkconfigure(entry.Name2,state = "disabled")
				tkconfigure(cb3,state = "normal")
				tkconfigure(entry.Name3,state = "normal")
				tkconfigure(lab30,state = "normal")
				tkconfigure(rb1,state = "normal")
				tkconfigure(rb2,state = "normal")
				tkconfigure(rb3,state = "normal")
				tkconfigure(rb4,state = "normal")
				tkconfigure(rb5,state = "normal")
				tkconfigure(aa2,state = "disabled")
				tkconfigure(cb5,state = "disabled")
				tkconfigure(cb6,state = "disabled")
				tkconfigure(entry.Name6,state = "disabled")
				tkconfigure(cb7,state = "disabled")
				tkconfigure(cb7.2,state = "disabled")
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
				tkconfigure(lab7.2,state = "normal")
				tkconfigure(lab8,state = "normal")
				tkconfigure(lab8.1,state = "normal")
				if(tclvalue(cbValue9)=="0" && tclvalue(cbValue10)=="0" && tclvalue(cbValue11)=="0" && tclvalue(cbValue12)=="0"){
					tkconfigure(entry.Name8,state = "disabled")
				}else{
					tkconfigure(entry.Name8,state = "normal")
				}
				tkconfigure(cb9,state = "normal")
				tkconfigure(lab9.1,state = "normal")
				tkconfigure(cb10,state = "normal")
				tkconfigure(cb11,state = "normal")
				tkconfigure(cb12,state = "normal")
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
				tkconfigure(OK.but,state = "normal")
			}
			else if( (tclvalue(cbValue1)=="1" && tclvalue(cbValue2)=="0" && tclvalue(cbValue3)=="0" ) || 
			    (tclvalue(cbValue1)=="0" && tclvalue(cbValue2)=="0" && tclvalue(cbValue3)=="0" )){
				tkconfigure(lab01,state = "normal")
				tkconfigure(lab02,state = "normal")
				tkconfigure(cb0,state = "normal")
				tkconfigure(entry.Name0,state = "normal")
				tkconfigure(lab00,state = "normal")
				tkconfigure(cb1,state = "normal")
				if(tclvalue(cbValue1)=="1"){
					tkconfigure(entry.Name1,state = "normal")
				}else{
					tkconfigure(entry.Name1,state = "disabled")
				}
				tkconfigure(cb2,state = "normal")
				if(tclvalue(cbValue2)=="1"){
					tkconfigure(entry.Name2,state = "normal")
				}else{
					tkconfigure(entry.Name2,state = "disabled")
				}
				tkconfigure(cb3,state = "normal")
				if(tclvalue(cbValue3)=="1"){
					tkconfigure(entry.Name3,state = "normal")
				}else{
					tkconfigure(entry.Name3,state = "disabled")
				}
				tkconfigure(lab30,state = "disabled")
				tkconfigure(rb1,state = "disabled")
				tkconfigure(rb2,state = "disabled")
				tkconfigure(rb3,state = "disabled")
				tkconfigure(rb4,state = "disabled")
				tkconfigure(rb5,state = "disabled")
				tkconfigure(aa2,state = "disabled")
				tkconfigure(cb5,state = "disabled")
				tkconfigure(cb6,state = "disabled")
				tkconfigure(entry.Name6,state = "disabled")
				tkconfigure(cb7,state = "disabled")
				tkconfigure(cb7.2,state = "disabled")
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
				tkconfigure(lab7.2,state = "disabled")
				tkconfigure(lab8,state = "disabled")
				tkconfigure(lab8.1,state = "disabled")
				tkconfigure(entry.Name8,state = "disabled")
				tkconfigure(cb9,state = "disabled")
				tkconfigure(lab9.1,state = "disabled")
				tkconfigure(cb10,state = "disabled")
				tkconfigure(cb11,state = "disabled")
				tkconfigure(cb12,state = "disabled")
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
				tkconfigure(OK.but,state = "normal")
			}
			else{
				tkconfigure(lab01,state = "normal")
				tkconfigure(lab02,state = "normal")
				tkconfigure(cb0,state = "normal")
				tkconfigure(entry.Name0,state = "normal")
				tkconfigure(lab00,state = "normal")
				tkconfigure(cb1,state = "normal")
				if(tclvalue(cbValue1)=="1"){
					tkconfigure(entry.Name1,state = "normal")
				}else{
					tkconfigure(entry.Name1,state = "disabled")
				}
				tkconfigure(cb2,state = "normal")
				if(tclvalue(cbValue2)=="1"){
					tkconfigure(entry.Name2,state = "normal")
				}else{
					tkconfigure(entry.Name2,state = "disabled")
				}
				tkconfigure(cb3,state = "normal")
				if(tclvalue(cbValue3)=="1"){
					tkconfigure(entry.Name3,state = "normal")
				}else{
					tkconfigure(entry.Name3,state = "disabled")
				}
				tkconfigure(lab30,state = "normal")
				tkconfigure(rb1,state = "normal")
				tkconfigure(rb2,state = "normal")
				tkconfigure(rb3,state = "normal")
				tkconfigure(rb4,state = "normal")
				tkconfigure(rb5,state = "normal")
				tkconfigure(aa2,state = "normal")
				tkconfigure(cb5,state = "normal")
				tkconfigure(cb6,state = "normal")
				tkconfigure(cb7,state = "normal")
				tkconfigure(cb7.2,state = "normal")
				if(tclvalue(cbValue6)=="1"){
					tkconfigure(entry.Name6,state = "normal")
					tkconfigure(lab6.91,state = "normal")
					tkconfigure(lab6.92,state = "normal")
					tkconfigure(lab7.from,state = "normal")
					tkconfigure(lab7.to,state = "normal")
					tkconfigure(lab7.incr,state = "normal")
					tkconfigure(entry.Name7.1,state = "normal")
					tkconfigure(entry.Name7.2,state = "normal")
					tkconfigure(entry.Name7.3,state = "normal")
				}
				if(tclvalue(cbValue7)=="1"){
					tkconfigure(lab6.91,state = "normal")
					tkconfigure(lab6.92,state = "normal")
					tkconfigure(lab7.from,state = "normal")
					tkconfigure(lab7.to,state = "normal")
					tkconfigure(lab7.incr,state = "normal")
					tkconfigure(entry.Name7.1,state = "normal")
					tkconfigure(entry.Name7.2,state = "normal")
					tkconfigure(entry.Name7.3,state = "normal")
					tkconfigure(lab13,state = "normal")
					tkconfigure(cb14,state = "normal")
					if(tclvalue(cbValue14)=="1") {
						tkconfigure(lab14,state = "normal")
						tkconfigure(cb15_1,state = "normal")
						tkconfigure(cb15_2,state = "normal")
						tkconfigure(lab15,state = "normal")
						tkconfigure(cb15,state = "normal")
						tkconfigure(cb16,state = "normal")
					}
				}
				if(tclvalue(cbValue7.2)=="1"){
					tkconfigure(lab6.91,state = "normal")
					tkconfigure(lab6.92,state = "normal")
					tkconfigure(lab7.from,state = "normal")
					tkconfigure(lab7.to,state = "normal")
					tkconfigure(lab7.incr,state = "normal")
					tkconfigure(entry.Name7.1,state = "normal")
					tkconfigure(entry.Name7.2,state = "normal")
					tkconfigure(entry.Name7.3,state = "normal")
					tkconfigure(lab13,state = "normal")
					tkconfigure(cb14,state = "normal")
					if(tclvalue(cbValue14)=="1") {
						tkconfigure(lab14,state = "normal")
						tkconfigure(cb15_1,state = "normal")
						tkconfigure(cb15_2,state = "normal")
						tkconfigure(lab15,state = "normal")
						tkconfigure(cb15,state = "normal")
						tkconfigure(cb16,state = "normal")
					}
				}
				tkconfigure(lab7.2,state = "normal")
				tkconfigure(lab8,state = "normal")
				tkconfigure(lab8.1,state = "normal")
				if(tclvalue(cbValue9)=="0" && tclvalue(cbValue10)=="0" && tclvalue(cbValue11)=="0" && tclvalue(cbValue12)=="0"){
					tkconfigure(entry.Name8,state = "disabled")
				}else{
					tkconfigure(entry.Name8,state = "normal")
				}
				tkconfigure(cb9,state = "normal")
				tkconfigure(lab9.1,state = "normal")
				tkconfigure(cb10,state = "normal")
				tkconfigure(cb11,state = "normal")
				tkconfigure(cb12,state = "normal")
				tkconfigure(OK.but,state = "normal")
			}
	    }
	    if(tclvalue(cbValue00)=="1"){
			tkconfigure(sdma,state = "disabled")
			tkconfigure(sdma2,state = "disabled")
			tkconfigure(lab01,state = "disabled")
			tkconfigure(lab02,state = "disabled")
			tkconfigure(cb0,state = "disabled")
			tkconfigure(entry.Name0,state = "disabled")
			tkconfigure(lab00,state = "disabled")
			tkconfigure(cb1,state = "disabled")
			tkconfigure(entry.Name1,state = "disabled")
			tkconfigure(cb2,state = "disabled")
			tkconfigure(entry.Name2,state = "disabled")
			tkconfigure(cb3,state = "disabled")
			tkconfigure(entry.Name3,state = "disabled")
			tkconfigure(lab30,state = "disabled")
			tkconfigure(rb1,state = "disabled")
			tkconfigure(rb2,state = "disabled")
			tkconfigure(rb3,state = "disabled")
			tkconfigure(rb4,state = "disabled")
			tkconfigure(rb5,state = "disabled")
			tkconfigure(aa2,state = "disabled")
			tkconfigure(cb5,state = "disabled")
			tkconfigure(cb6,state = "disabled")
			tkconfigure(entry.Name6,state = "disabled")
			tkconfigure(cb7,state = "disabled")
			tkconfigure(cb7.2,state = "disabled")
			tkconfigure(lab6.91,state = "disabled")
			tkconfigure(lab6.92,state = "disabled")
			tkconfigure(lab7.from,state = "disabled")
			tkconfigure(lab7.to,state = "disabled")
			tkconfigure(lab7.incr,state = "disabled")
			tkconfigure(entry.Name7.1,state = "disabled")
			tkconfigure(entry.Name7.2,state = "disabled")
			tkconfigure(entry.Name7.3,state = "disabled")
			tkconfigure(lab7.2,state = "disabled")
			tkconfigure(lab8,state = "disabled")
			tkconfigure(lab8.1,state = "disabled")
			tkconfigure(entry.Name8,state = "disabled")
			tkconfigure(cb9,state = "disabled")
			tkconfigure(lab9.1,state = "disabled")
			tkconfigure(cb10,state = "disabled")
			tkconfigure(cb11,state = "disabled")
			tkconfigure(cb12,state = "disabled")
			tkconfigure(lab13,state = "disabled")
			tkconfigure(cb14,state = "disabled")
			tkconfigure(lab14,state = "disabled")
			tkconfigure(cb15_1,state = "disabled")
			tkconfigure(cb15_2,state = "disabled")
			tkconfigure(lab15,state = "disabled")
			tkconfigure(cb15,state = "disabled")
			tkconfigure(cb16,state = "disabled")
			tkconfigure(OK.but,state = "disabled")
	    }
    }
    tkbind(wiz.cb,"<1>",set.entry.state.wiz)

	set.entry.state.seta<-function(){
		if(tclvalue(cbValue1)=="1"){
			if(tclvalue(cbValue00)=="1") tkconfigure(entry.Name1,state = "disabled")
			if(tclvalue(cbValue2)=="0"){
				tkconfigure(aa2,state = "disabled")
				tkconfigure(cb5,state = "disabled")
				tkconfigure(cb6,state = "disabled")
				tkconfigure(entry.Name6,state = "disabled")
				tkconfigure(cb7,state = "disabled")
				tkconfigure(cb7.2,state = "disabled")
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
			}
			if(tclvalue(cbValue2)=="0" && tclvalue(cbValue3)=="0"){
				tkconfigure(lab30,state = "disabled")
				tkconfigure(rb1,state = "disabled")
				tkconfigure(rb2,state = "disabled")
				tkconfigure(rb3,state = "disabled")
				tkconfigure(rb4,state = "disabled")
				tkconfigure(rb5,state = "disabled")
				tkconfigure(aa2,state = "disabled")
				tkconfigure(cb5,state = "disabled")
				tkconfigure(cb6,state = "disabled")
				tkconfigure(entry.Name6,state = "disabled")
				tkconfigure(cb7,state = "disabled")
				tkconfigure(cb7.2,state = "disabled")
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
				tkconfigure(lab7.2,state = "disabled")
				tkconfigure(lab8,state = "disabled")
				tkconfigure(lab8.1,state = "disabled")
				tkconfigure(entry.Name8,state = "disabled")
				tkconfigure(cb9,state = "disabled")
				tkconfigure(lab9.1,state = "disabled")
				tkconfigure(cb10,state = "disabled")
				tkconfigure(cb11,state = "disabled")
				tkconfigure(cb12,state = "disabled")
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
			}
		}
		if(tclvalue(cbValue1)=="0" && tclvalue(cbValue00)=="1"){
			if(tclvalue(cbValue2)=="0" && tclvalue(cbValue3)=="0"){
				tkconfigure(entry.Name1,state = "normal")
				tkconfigure(lab30,state = "disabled")
				tkconfigure(rb1,state = "disabled")
				tkconfigure(rb2,state = "disabled")
				tkconfigure(rb3,state = "disabled")
				tkconfigure(rb4,state = "disabled")
				tkconfigure(rb5,state = "disabled")
				tkconfigure(aa2,state = "disabled")
				tkconfigure(cb5,state = "disabled")
				tkconfigure(cb6,state = "disabled")
				tkconfigure(entry.Name6,state = "disabled")
				tkconfigure(cb7,state = "disabled")
				tkconfigure(cb7.2,state = "disabled")
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
				tkconfigure(lab7.2,state = "disabled")
				tkconfigure(lab8,state = "disabled")
				tkconfigure(lab8.1,state = "disabled")
				tkconfigure(entry.Name8,state = "disabled")
				tkconfigure(cb9,state = "disabled")
				tkconfigure(lab9.1,state = "disabled")
				tkconfigure(cb10,state = "disabled")
				tkconfigure(cb11,state = "disabled")
				tkconfigure(cb12,state = "disabled")
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
			}else{
				tkconfigure(entry.Name1,state = "normal")
				tkconfigure(lab30,state = "normal")
				tkconfigure(rb1,state = "normal")
				tkconfigure(rb2,state = "normal")
				tkconfigure(rb3,state = "normal")
				tkconfigure(rb4,state = "normal")
				tkconfigure(rb5,state = "normal")
				tkconfigure(aa2,state = "normal")
				tkconfigure(cb5,state = "normal")
				tkconfigure(cb6,state = "normal")
				tkconfigure(cb7,state = "normal")
				tkconfigure(cb7.2,state = "normal")
				if(tclvalue(cbValue6)=="1"){
					tkconfigure(entry.Name6,state = "normal")
					tkconfigure(lab6.91,state = "normal")
					tkconfigure(lab6.92,state = "normal")
					tkconfigure(lab7.from,state = "normal")
					tkconfigure(lab7.to,state = "normal")
					tkconfigure(lab7.incr,state = "normal")
					tkconfigure(entry.Name7.1,state = "normal")
					tkconfigure(entry.Name7.2,state = "normal")
					tkconfigure(entry.Name7.3,state = "normal")
				}
				if(tclvalue(cbValue7)=="1"){
					tkconfigure(lab6.91,state = "normal")
					tkconfigure(lab6.92,state = "normal")
					tkconfigure(lab7.from,state = "normal")
					tkconfigure(lab7.to,state = "normal")
					tkconfigure(lab7.incr,state = "normal")
					tkconfigure(entry.Name7.1,state = "normal")
					tkconfigure(entry.Name7.2,state = "normal")
					tkconfigure(entry.Name7.3,state = "normal")
					tkconfigure(lab13,state = "normal")
					tkconfigure(cb14,state = "normal")
					if(tclvalue(cbValue14)=="1") {
						tkconfigure(lab14,state = "normal")
						tkconfigure(cb15_1,state = "normal")
						tkconfigure(cb15_2,state = "normal")
						tkconfigure(lab15,state = "normal")
						tkconfigure(cb15,state = "normal")
						tkconfigure(cb16,state = "normal")
					}
				}
				if(tclvalue(cbValue7.2)=="1"){
					tkconfigure(lab6.91,state = "normal")
					tkconfigure(lab6.92,state = "normal")
					tkconfigure(lab7.from,state = "normal")
					tkconfigure(lab7.to,state = "normal")
					tkconfigure(lab7.incr,state = "normal")
					tkconfigure(entry.Name7.1,state = "normal")
					tkconfigure(entry.Name7.2,state = "normal")
					tkconfigure(entry.Name7.3,state = "normal")
					tkconfigure(lab13,state = "normal")
					tkconfigure(cb14,state = "normal")
					if(tclvalue(cbValue14)=="1") {
						tkconfigure(lab14,state = "normal")
						tkconfigure(cb15_1,state = "normal")
						tkconfigure(cb15_2,state = "normal")
						tkconfigure(lab15,state = "normal")
						tkconfigure(cb15,state = "normal")
						tkconfigure(cb16,state = "normal")
					}
				}
			}
		}
	}
	tkbind(cb1,"<1>",set.entry.state.seta)

	set.entry.state.sb<-function(){
		if(tclvalue(cbValue2)=="1"){
			if(tclvalue(cbValue00)=="1") tkconfigure(entry.Name2,state = "disabled")
			if(tclvalue(cbValue1)=="0" && tclvalue(cbValue3)=="1"){
				tkconfigure(aa2,state = "disabled")
				tkconfigure(cb5,state = "disabled")
				tkconfigure(cb6,state = "disabled")
				tkconfigure(entry.Name6,state = "disabled")
				tkconfigure(cb7,state = "disabled")
				tkconfigure(cb7.2,state = "disabled")
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
			}
			if(tclvalue(cbValue1)=="0" && tclvalue(cbValue3)=="0"){
				tkconfigure(lab30,state = "disabled")
				tkconfigure(rb1,state = "disabled")
				tkconfigure(rb2,state = "disabled")
				tkconfigure(rb3,state = "disabled")
				tkconfigure(rb4,state = "disabled")
				tkconfigure(rb5,state = "disabled")
				tkconfigure(aa2,state = "disabled")
				tkconfigure(cb5,state = "disabled")
				tkconfigure(cb6,state = "disabled")
				tkconfigure(entry.Name6,state = "disabled")
				tkconfigure(cb7,state = "disabled")
				tkconfigure(cb7.2,state = "disabled")
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
				tkconfigure(lab7.2,state = "disabled")
				tkconfigure(lab8,state = "disabled")
				tkconfigure(lab8.1,state = "disabled")
				tkconfigure(entry.Name8,state = "disabled")
				tkconfigure(cb9,state = "disabled")
				tkconfigure(lab9.1,state = "disabled")
				tkconfigure(cb10,state = "disabled")
				tkconfigure(cb11,state = "disabled")
				tkconfigure(cb12,state = "disabled")
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
			}
			if(tclvalue(cbValue1)=="1" && tclvalue(cbValue3)=="0"){
				tkconfigure(lab30,state = "disabled")
				tkconfigure(rb1,state = "disabled")
				tkconfigure(rb2,state = "disabled")
				tkconfigure(rb3,state = "disabled")
				tkconfigure(rb4,state = "disabled")
				tkconfigure(rb5,state = "disabled")
				tkconfigure(aa2,state = "disabled")
				tkconfigure(cb5,state = "disabled")
				tkconfigure(cb6,state = "disabled")
				tkconfigure(entry.Name6,state = "disabled")
				tkconfigure(cb7,state = "disabled")
				tkconfigure(cb7.2,state = "disabled")
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
				tkconfigure(lab7.2,state = "disabled")
				tkconfigure(lab8,state = "disabled")
				tkconfigure(lab8.1,state = "disabled")
				tkconfigure(entry.Name8,state = "disabled")
				tkconfigure(cb9,state = "disabled")
				tkconfigure(lab9.1,state = "disabled")
				tkconfigure(cb10,state = "disabled")
				tkconfigure(cb11,state = "disabled")
				tkconfigure(cb12,state = "disabled")
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
			}
		}
		if(tclvalue(cbValue2)=="0" && tclvalue(cbValue00)=="1"){
			tkconfigure(lab30,state = "normal")
			tkconfigure(rb1,state = "normal")
			tkconfigure(rb2,state = "normal")
			tkconfigure(rb3,state = "normal")
			tkconfigure(rb4,state = "normal")
			tkconfigure(rb5,state = "normal")
			tkconfigure(entry.Name2,state = "normal")
			tkconfigure(aa2,state = "normal")
			tkconfigure(cb5,state = "normal")
			tkconfigure(cb6,state = "normal")
			tkconfigure(cb7,state = "normal")
			tkconfigure(cb7.2,state = "normal")
			if(tclvalue(cbValue6)=="1"){
				tkconfigure(entry.Name6,state = "normal")
				tkconfigure(lab6.91,state = "normal")
				tkconfigure(lab6.92,state = "normal")
				tkconfigure(lab7.from,state = "normal")
				tkconfigure(lab7.to,state = "normal")
				tkconfigure(lab7.incr,state = "normal")
				tkconfigure(entry.Name7.1,state = "normal")
				tkconfigure(entry.Name7.2,state = "normal")
				tkconfigure(entry.Name7.3,state = "normal")
			}
			if(tclvalue(cbValue7)=="1"){
				tkconfigure(lab6.91,state = "normal")
				tkconfigure(lab6.92,state = "normal")
				tkconfigure(lab7.from,state = "normal")
				tkconfigure(lab7.to,state = "normal")
				tkconfigure(lab7.incr,state = "normal")
				tkconfigure(entry.Name7.1,state = "normal")
				tkconfigure(entry.Name7.2,state = "normal")
				tkconfigure(entry.Name7.3,state = "normal")
				tkconfigure(lab13,state = "normal")
				tkconfigure(cb14,state = "normal")
				if(tclvalue(cbValue14)=="1") {
					tkconfigure(lab14,state = "normal")
					tkconfigure(cb15_1,state = "normal")
					tkconfigure(cb15_2,state = "normal")
					tkconfigure(lab15,state = "normal")
					tkconfigure(cb15,state = "normal")
					tkconfigure(cb16,state = "normal")
				}
			}
			tkconfigure(lab7.2,state = "normal")
			if(tclvalue(cbValue7.2)=="1"){
				tkconfigure(lab6.91,state = "normal")
				tkconfigure(lab6.92,state = "normal")
				tkconfigure(lab7.from,state = "normal")
				tkconfigure(lab7.to,state = "normal")
				tkconfigure(lab7.incr,state = "normal")
				tkconfigure(entry.Name7.1,state = "normal")
				tkconfigure(entry.Name7.2,state = "normal")
				tkconfigure(entry.Name7.3,state = "normal")
				tkconfigure(lab13,state = "normal")
				tkconfigure(cb14,state = "normal")
				if(tclvalue(cbValue14)=="1") {
					tkconfigure(lab14,state = "normal")
					tkconfigure(cb15_1,state = "normal")
					tkconfigure(cb15_2,state = "normal")
					tkconfigure(lab15,state = "normal")
					tkconfigure(cb15,state = "normal")
					tkconfigure(cb16,state = "normal")
				}
			}
			tkconfigure(lab8,state = "normal")
			tkconfigure(lab8.1,state = "normal")
			if(tclvalue(cbValue9)=="0" && tclvalue(cbValue10)=="0" && tclvalue(cbValue11)=="0" && tclvalue(cbValue12)=="0"){
				tkconfigure(entry.Name8,state = "disabled")
			}else{
				tkconfigure(entry.Name8,state = "normal")
			}
			tkconfigure(cb9,state = "normal")
			tkconfigure(lab9.1,state = "normal")
			tkconfigure(cb10,state = "normal")
			tkconfigure(cb11,state = "normal")
			tkconfigure(cb12,state = "normal")
		}
	}
	tkbind(cb2,"<1>",set.entry.state.sb)

	set.entry.state.se<-function(){
		if(tclvalue(cbValue3)=="1"){
			if(tclvalue(cbValue00)=="1") tkconfigure(entry.Name3,state = "disabled")
			if(tclvalue(cbValue1)=="0" && tclvalue(cbValue2)=="0"){
				tkconfigure(lab30,state = "disabled")
				tkconfigure(rb1,state = "disabled")
				tkconfigure(rb2,state = "disabled")
				tkconfigure(rb3,state = "disabled")
				tkconfigure(rb4,state = "disabled")
				tkconfigure(rb5,state = "disabled")
				tkconfigure(lab7.2,state = "disabled")
				tkconfigure(lab8,state = "disabled")
				tkconfigure(lab8.1,state = "disabled")
				tkconfigure(entry.Name8,state = "disabled")
				tkconfigure(cb9,state = "disabled")
				tkconfigure(lab9.1,state = "disabled")
				tkconfigure(cb10,state = "disabled")
				tkconfigure(cb11,state = "disabled")
				tkconfigure(cb12,state = "disabled")
			}
			if(tclvalue(cbValue1)=="1" && tclvalue(cbValue2)=="0"){
				tkconfigure(lab30,state = "disabled")
				tkconfigure(rb1,state = "disabled")
				tkconfigure(rb2,state = "disabled")
				tkconfigure(rb3,state = "disabled")
				tkconfigure(rb4,state = "disabled")
				tkconfigure(rb5,state = "disabled")
				tkconfigure(aa2,state = "disabled")
				tkconfigure(cb5,state = "disabled")
				tkconfigure(cb6,state = "disabled")
				tkconfigure(entry.Name6,state = "disabled")
				tkconfigure(cb7,state = "disabled")
				tkconfigure(cb7.2,state = "disabled")
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
				tkconfigure(lab7.2,state = "disabled")
				tkconfigure(lab8,state = "disabled")
				tkconfigure(lab8.1,state = "disabled")
				tkconfigure(entry.Name8,state = "disabled")
				tkconfigure(cb9,state = "disabled")
				tkconfigure(lab9.1,state = "disabled")
				tkconfigure(cb10,state = "disabled")
				tkconfigure(cb11,state = "disabled")
				tkconfigure(cb12,state = "disabled")
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
			}
		}
		if(tclvalue(cbValue3)=="0" && tclvalue(cbValue00)=="1") {
			if(tclvalue(cbValue1)=="0" && tclvalue(cbValue2)=="0"){
				tkconfigure(entry.Name3,state = "normal")
				tkconfigure(lab30,state = "normal")
				tkconfigure(rb1,state = "normal")
				tkconfigure(rb2,state = "normal")
				tkconfigure(rb3,state = "normal")
				tkconfigure(rb4,state = "normal")
				tkconfigure(rb5,state = "normal")
				tkconfigure(lab7.2,state = "normal")
				tkconfigure(lab8,state = "normal")
				tkconfigure(lab8.1,state = "normal")
				if(tclvalue(cbValue9)=="0" && tclvalue(cbValue10)=="0" && tclvalue(cbValue11)=="0" && tclvalue(cbValue12)=="0"){
					tkconfigure(entry.Name8,state = "disabled")
				}else{
					tkconfigure(entry.Name8,state = "normal")
				}
				tkconfigure(cb9,state = "normal")
				tkconfigure(lab9.1,state = "normal")
				tkconfigure(cb10,state = "normal")
				tkconfigure(cb11,state = "normal")
				tkconfigure(cb12,state = "normal")
			}else{
				tkconfigure(entry.Name3,state = "normal")
				tkconfigure(lab30,state = "normal")
				tkconfigure(rb1,state = "normal")
				tkconfigure(rb2,state = "normal")
				tkconfigure(rb3,state = "normal")
				tkconfigure(rb4,state = "normal")
				tkconfigure(rb5,state = "normal")
				tkconfigure(aa2,state = "normal")
				tkconfigure(cb5,state = "normal")
				tkconfigure(cb6,state = "normal")
				tkconfigure(cb7,state = "normal")
				tkconfigure(cb7.2,state = "normal")
				if(tclvalue(cbValue6)=="1"){
					tkconfigure(entry.Name6,state = "normal")
					tkconfigure(lab6.91,state = "normal")
					tkconfigure(lab6.92,state = "normal")
					tkconfigure(lab7.from,state = "normal")
					tkconfigure(lab7.to,state = "normal")
					tkconfigure(lab7.incr,state = "normal")
					tkconfigure(entry.Name7.1,state = "normal")
					tkconfigure(entry.Name7.2,state = "normal")
					tkconfigure(entry.Name7.3,state = "normal")
				}
				if(tclvalue(cbValue7)=="1"){
					tkconfigure(lab6.91,state = "normal")
					tkconfigure(lab6.92,state = "normal")
					tkconfigure(lab7.from,state = "normal")
					tkconfigure(lab7.to,state = "normal")
					tkconfigure(lab7.incr,state = "normal")
					tkconfigure(entry.Name7.1,state = "normal")
					tkconfigure(entry.Name7.2,state = "normal")
					tkconfigure(entry.Name7.3,state = "normal")
					tkconfigure(lab13,state = "normal")
					tkconfigure(cb14,state = "normal")
					if(tclvalue(cbValue14)=="1") {
						tkconfigure(lab14,state = "normal")
						tkconfigure(cb15_1,state = "normal")
						tkconfigure(cb15_2,state = "normal")
						tkconfigure(lab15,state = "normal")
						tkconfigure(cb15,state = "normal")
						tkconfigure(cb16,state = "normal")
					}
				}
				if(tclvalue(cbValue7.2)=="1"){
					tkconfigure(lab6.91,state = "normal")
					tkconfigure(lab6.92,state = "normal")
					tkconfigure(lab7.from,state = "normal")
					tkconfigure(lab7.to,state = "normal")
					tkconfigure(lab7.incr,state = "normal")
					tkconfigure(entry.Name7.1,state = "normal")
					tkconfigure(entry.Name7.2,state = "normal")
					tkconfigure(entry.Name7.3,state = "normal")
					tkconfigure(lab13,state = "normal")
					tkconfigure(cb14,state = "normal")
					if(tclvalue(cbValue14)=="1") {
						tkconfigure(lab14,state = "normal")
						tkconfigure(cb15_1,state = "normal")
						tkconfigure(cb15_2,state = "normal")
						tkconfigure(lab15,state = "normal")
						tkconfigure(cb15,state = "normal")
						tkconfigure(cb16,state = "normal")
					}
				}
				tkconfigure(lab7.2,state = "normal")
				tkconfigure(lab8,state = "normal")
				tkconfigure(lab8.1,state = "normal")
				if(tclvalue(cbValue9)=="0" && tclvalue(cbValue10)=="0" && tclvalue(cbValue11)=="0" && tclvalue(cbValue12)=="0"){
					tkconfigure(entry.Name8,state = "disabled")
				}else{
					tkconfigure(entry.Name8,state = "normal")
				}
				tkconfigure(cb9,state = "normal")
				tkconfigure(lab9.1,state = "normal")
				tkconfigure(cb10,state = "normal")
				tkconfigure(cb11,state = "normal")
				tkconfigure(cb12,state = "normal")
			}
		}
	}
	tkbind(cb3,"<1>",set.entry.state.se)

	set.entry.state.q<-function(){
		if(tclvalue(cbValue6)=="1") {
			tkconfigure(entry.Name6,state = "disabled")
			if(tclvalue(cbValue7)=="0" && tclvalue(cbValue7.2)=="0"){
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
					tkconfigure(lab7.from,state = "disabled")
					tkconfigure(lab7.to,state = "disabled")
					tkconfigure(lab7.incr,state = "disabled")
					tkconfigure(entry.Name7.1,state = "disabled")
					tkconfigure(entry.Name7.2,state = "disabled")
					tkconfigure(entry.Name7.3,state = "disabled")
			}
		}
		if(tclvalue(cbValue6)=="0") {
			if( ( tclvalue(cbValue2)=="1" || ( tclvalue(cbValue1)=="1" && tclvalue(cbValue3)=="1" ) ) && tclvalue(cbValue00)=="1" ) {
				tkconfigure(entry.Name6,state = "normal")
				tkconfigure(lab6.91,state = "normal")
				tkconfigure(lab6.92,state = "normal")
				tkconfigure(lab7.from,state = "normal")
				tkconfigure(lab7.to,state = "normal")
				tkconfigure(lab7.incr,state = "normal")
				tkconfigure(entry.Name7.1,state = "normal")
				tkconfigure(entry.Name7.2,state = "normal")
				tkconfigure(entry.Name7.3,state = "normal")
			}
		}
	}
	tkbind(cb6,"<1>",set.entry.state.q)

	set.entry.state.cdf<-function(){
		if(tclvalue(cbValue7)=="1") {
			if(tclvalue(cbValue7.2)=="0"){
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
			}
			if(tclvalue(cbValue6)=="0" && tclvalue(cbValue7.2)=="0"){
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
			}

		}
		if(tclvalue(cbValue7)=="0") {
			if( ( tclvalue(cbValue2)=="1" || ( tclvalue(cbValue1)=="1" && tclvalue(cbValue3)=="1" ) ) && tclvalue(cbValue00)=="1" ) {
				tkconfigure(lab6.91,state = "normal")
				tkconfigure(lab6.92,state = "normal")
				tkconfigure(lab7.from,state = "normal")
				tkconfigure(lab7.to,state = "normal")
				tkconfigure(lab7.incr,state = "normal")
				tkconfigure(entry.Name7.1,state = "normal")
				tkconfigure(entry.Name7.2,state = "normal")
				tkconfigure(entry.Name7.3,state = "normal")
				tkconfigure(lab13,state = "normal")
				tkconfigure(cb14,state = "normal")
				if(tclvalue(cbValue14)=="1") {
					tkconfigure(lab14,state = "normal")
					tkconfigure(cb15_1,state = "normal")
					tkconfigure(cb15_2,state = "normal")
					tkconfigure(lab15,state = "normal")
					tkconfigure(cb15,state = "normal")
					tkconfigure(cb16,state = "normal")
				}
			}
		}
	}
	tkbind(cb7,"<1>",set.entry.state.cdf)

	set.entry.state.pdf<-function(){
		if(tclvalue(cbValue7.2)=="1") {
			if(tclvalue(cbValue7)=="0"){
				tkconfigure(lab13,state = "disabled")
				tkconfigure(cb14,state = "disabled")
				tkconfigure(lab14,state = "disabled")
				tkconfigure(cb15_1,state = "disabled")
				tkconfigure(cb15_2,state = "disabled")
				tkconfigure(lab15,state = "disabled")
				tkconfigure(cb15,state = "disabled")
				tkconfigure(cb16,state = "disabled")
			}
			if(tclvalue(cbValue6)=="0" && tclvalue(cbValue7)=="0"){
				tkconfigure(lab6.91,state = "disabled")
				tkconfigure(lab6.92,state = "disabled")
				tkconfigure(lab7.from,state = "disabled")
				tkconfigure(lab7.to,state = "disabled")
				tkconfigure(lab7.incr,state = "disabled")
				tkconfigure(entry.Name7.1,state = "disabled")
				tkconfigure(entry.Name7.2,state = "disabled")
				tkconfigure(entry.Name7.3,state = "disabled")
			}
		}
		if(tclvalue(cbValue7.2)=="0") {
			if( ( tclvalue(cbValue2)=="1" || ( tclvalue(cbValue1)=="1" && tclvalue(cbValue3)=="1" ) ) && tclvalue(cbValue00)=="1" ) {
				tkconfigure(lab6.91,state = "normal")
				tkconfigure(lab6.92,state = "normal")
				tkconfigure(lab7.from,state = "normal")
				tkconfigure(lab7.to,state = "normal")
				tkconfigure(lab7.incr,state = "normal")
				tkconfigure(entry.Name7.1,state = "normal")
				tkconfigure(entry.Name7.2,state = "normal")
				tkconfigure(entry.Name7.3,state = "normal")
				tkconfigure(lab13,state = "normal")
				tkconfigure(cb14,state = "normal")
				if(tclvalue(cbValue14)=="1") {
					tkconfigure(lab14,state = "normal")
					tkconfigure(cb15_1,state = "normal")
					tkconfigure(cb15_2,state = "normal")
					tkconfigure(lab15,state = "normal")
					tkconfigure(cb15,state = "normal")
					tkconfigure(cb16,state = "normal")
				}
			}
		}
	}
	tkbind(cb7.2,"<1>",set.entry.state.pdf)

	set.entry.state.fim<-function(){
		if(tclvalue(cbValue9)=="1") {
			if(tclvalue(cbValue10)=="0" && tclvalue(cbValue11)=="0" && tclvalue(cbValue12)=="0"){
				tkconfigure(entry.Name8,state = "disabled")
			}
		}
		if( tclvalue(cbValue9)=="0" && tclvalue(cbValue00)=="1" && ( tclvalue(cbValue2)=="1" || tclvalue(cbValue3)=="1" ) ) {
			tkconfigure(entry.Name8,state = "normal")
		}
	}
	tkbind(cb9,"<1>",set.entry.state.fim)


	set.entry.state.oima<-function(){
		if(tclvalue(cbValue10)=="1") {
			if(tclvalue(cbValue9)=="0" && tclvalue(cbValue11)=="0" && tclvalue(cbValue12)=="0"){
				tkconfigure(entry.Name8,state = "disabled")
			}
		}
		if( tclvalue(cbValue10)=="0" && tclvalue(cbValue00)=="1" && ( tclvalue(cbValue2)=="1" || tclvalue(cbValue3)=="1" ) ) {
			tkconfigure(entry.Name8,state = "normal")
		}
	}
	tkbind(cb10,"<1>",set.entry.state.oima)


	set.entry.state.oimb<-function(){
		if(tclvalue(cbValue11)=="1") {
			if(tclvalue(cbValue9)=="0" && tclvalue(cbValue10)=="0" && tclvalue(cbValue12)=="0"){
				tkconfigure(entry.Name8,state = "disabled")
			}
		}
		if( tclvalue(cbValue11)=="0" && tclvalue(cbValue00)=="1" && ( tclvalue(cbValue2)=="1" || tclvalue(cbValue3)=="1" ) ) {
			tkconfigure(entry.Name8,state = "normal")
		}
	}
	tkbind(cb11,"<1>",set.entry.state.oimb)


	set.entry.state.oimc<-function(){
		if(tclvalue(cbValue12)=="1") {
			if(tclvalue(cbValue9)=="0" && tclvalue(cbValue10)=="0" && tclvalue(cbValue11)=="0"){
				tkconfigure(entry.Name8,state = "disabled")
			}
		}
		if( tclvalue(cbValue12)=="0" && tclvalue(cbValue00)=="1" && ( tclvalue(cbValue2)=="1" || tclvalue(cbValue3)=="1" ) ) {
			tkconfigure(entry.Name8,state = "normal")
		}
	}
	tkbind(cb12,"<1>",set.entry.state.oimc)


	set.entry.state.pseudo<-function(){
		if(tclvalue(cbValue14)=="1") {
			tkconfigure(lab14,state = "disabled")
			tkconfigure(cb15_1,state = "disabled")
			tkconfigure(cb15_2,state = "disabled")
			tkconfigure(lab15,state = "disabled")
			tkconfigure(cb15,state = "disabled")
			tkconfigure(cb16,state = "disabled")
		}
		if(tclvalue(cbValue14)=="0") {
			if( (tclvalue(cbValue7)=="1" || tclvalue(cbValue7.2)=="1") && tclvalue(cbValue00)=="1" && ( tclvalue(cbValue2)=="1" || ( tclvalue(cbValue1)=="1" && tclvalue(cbValue3)=="1" ) ) ) {
				tkconfigure(lab14,state = "normal")
				tkconfigure(cb15_1,state = "normal")
				tkconfigure(cb15_2,state = "normal")
				tkconfigure(lab15,state = "normal")
				tkconfigure(cb15,state = "normal")
				tkconfigure(cb16,state = "normal")
			}
		}
	}
	tkbind(cb14,"<1>",set.entry.state.pseudo)

    pseudocheck <- tclVar("0")

    tkgrid(f1, sticky="w")
    tkgrid(f1.1, sticky="w")
    tkgrid(f1.2, sticky="w")
    tkgrid(f2, sticky="w")
    tkgrid(f3, sticky="w")
    tkgrid(f3.1, sticky="w")
    tkgrid(f3.2, sticky="w")


    OnOK <- function(){
	  
        if(!exists(tclvalue(datapick))){
		stop("No data set! Please choose one to analyze!")
		
        }
        wVal <- eval(parse(text=tclvalue(w)))
        qVal <- eval(parse(text=tclvalue(q)))
        alphaVal <- eval(parse(text=tclvalue(alpha)))
        etaVal <- eval(parse(text=tclvalue(eta.int)))
        setaVal <- eval(parse(text=tclvalue(seta.int)))
        sbVal <- eval(parse(text=tclvalue(sb.int)))
        seVal <- eval(parse(text=tclvalue(se.int)))
        pseudoVal <- as.character(tclvalue(cbValue14))
        cbVal_0 <- as.character(tclvalue(cbValue0))
        cbVal_1 <- as.character(tclvalue(cbValue1))
        cbVal_2 <- as.character(tclvalue(cbValue2))
        cbVal_3 <- as.character(tclvalue(cbValue3))
        MTTFVal <- as.character(tclvalue(cbValue5))
        quanVal <- as.character(tclvalue(cbValue6))
        FIMVal <- as.numeric(tclvalue(cbValue9))
        OIMAVal <- as.numeric(tclvalue(cbValue10))
        OIMBVal <- as.numeric(tclvalue(cbValue11))
        OIMCVal <- as.numeric(tclvalue(cbValue12))
        cdftVal <- seq( eval(parse(text=tclvalue(cdft1))) ,eval(parse(text=tclvalue(cdft2))), eval(parse(text=tclvalue(cdft3))) )
        hVal <- eval(parse(text=tclvalue(h)))
        pdf_plot <- as.character(tclvalue(cbValue7.2))
        cdf_plot <- as.character(tclvalue(cbValue7))
        ppplotVal <- as.character(tclvalue(cbValue15_1))
        qqplotVal <- as.character(tclvalue(cbValue15_2))
        kstestVal <- as.character(tclvalue(cbValue15))
        adtestVal <- as.character(tclvalue(cbValue16))
        medthodChoice <- as.character(tclvalue(rbValue))

        cat("Data:",tclvalue(datapick),"\n\n")
        cat("Please wait, analyzing degradation data...","\n\n")

        if (cbVal_0=="1" && cbVal_1=="0" && cbVal_2=="0" && cbVal_3=="0")
        cat("Warning message: please choose at least one variation effect in the degradation model!","\n\n")

	  if(cbVal_0==1){
		  if(!is.numeric(etaVal) || length(etaVal)>1) stop("The inputted initial value must be numeric and a scale, not a vector!")
	  }
	  if(cbVal_1==1){
		  if(!is.numeric(setaVal) || length(setaVal)>1 || setaVal<0 ) stop("The inputted initial value must be numeric and a positive value, not a vector!")
	  }
	  if(cbVal_2==1){
		  if(!is.numeric(sbVal) || length(sbVal)>1 || sbVal<0 ) stop("The inputted initial value must be numeric and a positive value, not a vector!")
	  }
	  if(cbVal_3==1){
		  if(!is.numeric(seVal) || length(seVal)>1 || seVal<0 ) stop("The inputted initial value must be numeric and a positive value, not a vector!")
	  }
	  if(!is.numeric(wVal) || length(wVal)>1 || wVal<=0 ) stop("The inputted threshold must be numeric and a positive value, not a vector!")
	  if( FIMVal!=0 || OIMAVal!=0 || OIMBVal!=0 || OIMCVal!=0 ){
			if(!is.numeric(alphaVal) || length(alphaVal)>1 || alphaVal>=1 || alphaVal<=0) stop("The inputted significant level must be numeric and a vector between 0 and 1!")
	  }
	  if(quanVal==1){
		  if(!is.numeric(qVal) || any(qVal<0) || any(qVal>1)) stop("The inputted qth quantile must be numeric and a vector between 0 and 1!")
	  }
	  if(!is.numeric(cdftVal) || any(cdftVal<0) ) stop("The inputted interval (domain) must be numeric and should not be negative!")




        if (cbVal_0=="1" && cbVal_1=="1" && cbVal_2=="1" && cbVal_3=="1"){
	        Model_M0(get(tclvalue(datapick)),etaVal,setaVal,sbVal,seVal,wVal,qVal,alphaVal,FIMVal,OIMAVal,OIMBVal,OIMCVal,cdftVal,hVal,pseudoVal,MTTFVal,quanVal,pdf_plot,cdf_plot,ppplotVal,qqplotVal,kstestVal,adtestVal,tclvalue(datapick),medthodChoice,e2)
		  tclvalue(eta.int) <- get('par.M0',envir=e2)[1]
		  tclvalue(seta.int) <- get('par.M0',envir=e2)[2]
		  tclvalue(sb.int) <- get('par.M0',envir=e2)[3]
		  tclvalue(se.int) <- get('par.M0',envir=e2)[4]
	  }
        if (cbVal_0=="1" && cbVal_1=="1" && cbVal_2=="0" && cbVal_3=="1"){
	        Model_M1(get(tclvalue(datapick)),etaVal,setaVal,seVal,wVal,qVal,alphaVal,FIMVal,OIMAVal,OIMBVal,OIMCVal,cdftVal,hVal,pseudoVal,MTTFVal,quanVal,pdf_plot,cdf_plot,ppplotVal,qqplotVal,kstestVal,adtestVal,tclvalue(datapick),medthodChoice,e2)
		  tclvalue(eta.int) <- get('par.M1',envir=e2)[1]
		  tclvalue(seta.int) <- get('par.M1',envir=e2)[2]
		  tclvalue(se.int) <- get('par.M1',envir=e2)[3]
	  }
        if (cbVal_0=="1" && cbVal_1=="0" && cbVal_2=="1" && cbVal_3=="0"){
	        Model_M2(get(tclvalue(datapick)),etaVal,sbVal,wVal,qVal,alphaVal,FIMVal,OIMAVal,OIMBVal,OIMCVal,cdftVal,pseudoVal,MTTFVal,quanVal,pdf_plot,cdf_plot,ppplotVal,qqplotVal,kstestVal,adtestVal,tclvalue(datapick),e2)
		  tclvalue(eta.int) <- get('par.M2',envir=e2)[1]
		  tclvalue(sb.int) <- get('par.M2',envir=e2)[2]
	  }
        if (cbVal_0=="1" && cbVal_1=="1" && cbVal_2=="1" && cbVal_3=="0"){
	        Model_M3(get(tclvalue(datapick)),etaVal,setaVal,sbVal,wVal,qVal,alphaVal,FIMVal,OIMAVal,OIMBVal,OIMCVal,cdftVal,hVal,pseudoVal,MTTFVal,quanVal,pdf_plot,cdf_plot,ppplotVal,qqplotVal,kstestVal,adtestVal,tclvalue(datapick),medthodChoice,e2)
		  tclvalue(eta.int) <- get('par.M3',envir=e2)[1]
		  tclvalue(seta.int) <- get('par.M3',envir=e2)[2]
		  tclvalue(sb.int) <- get('par.M3',envir=e2)[3]
	  }
        if (cbVal_0=="1" && cbVal_1=="0" && cbVal_2=="0" && cbVal_3=="1"){
	        Model_M4(get(tclvalue(datapick)),etaVal,seVal,wVal,qVal,alphaVal,FIMVal,OIMAVal,OIMBVal,OIMCVal,e2)
		  tclvalue(eta.int) <- get('par.M4',envir=e2)[1]
		  tclvalue(se.int) <- get('par.M4',envir=e2)[2]
	  }
        if (cbVal_0=="1" && cbVal_1=="0" && cbVal_2=="1" && cbVal_3=="1"){
	        Model_M5(get(tclvalue(datapick)),etaVal,sbVal,seVal,wVal,qVal,alphaVal,FIMVal,OIMAVal,OIMBVal,OIMCVal,cdftVal,pseudoVal,MTTFVal,quanVal,pdf_plot,cdf_plot,ppplotVal,qqplotVal,kstestVal,adtestVal,tclvalue(datapick),medthodChoice,e2)
		  tclvalue(eta.int) <- get('par.M5',envir=e2)[1]
		  tclvalue(sb.int) <- get('par.M5',envir=e2)[2]
		  tclvalue(se.int) <- get('par.M5',envir=e2)[3]
	  }
        if (cbVal_0=="1" && cbVal_1=="1" && cbVal_2=="0" && cbVal_3=="0")
        cat("Warning message: the covariance matrix of this degradation model is singular!","\n\n")

    }
    OnCancel <- function() {
        tkgrab.release(base)
        tkdestroy(base)
    }

    F5 <- tkframe(mainfrm,relief="groove")
    OK.but <- tkbutton(F5,text="Run",width=10,command=OnOK, state='disabled')
    Cancel.but <- tkbutton(F5,text="Exit",width=10,command=OnCancel)
    tkgrid(OK.but,tklabel(F5,text="		"),Cancel.but)
    tkgrid(Mainfrm, sticky="w")

    tkgrid(F)
    tkgrid(nb)
    tkgrid(F5)
    tkgrid(mainfrm)
    tkgrid(Mainfrm)

}

