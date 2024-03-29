lineshapes <- list(
	none=list(	
		func="(1i*(-1 + exp(-(1i*aq*(omega - omega0 - 1i*r2)))))/(omega - omega0 - 1i*r2)",
		domega0="(1i - 1i*exp(aq*(1i*(omega - omega0) + r2)) + aq*(-omega + omega0 + 1i*r2))/(exp(aq*(1i*(omega - omega0) + r2))*(-omega + omega0 + 1i*r2)^2)",
		dr2="(-1 + exp(aq*(1i*(omega - omega0) + r2)) - 1i*aq*(omega - omega0 - 1i*r2))/(exp(aq*(1i*(omega - omega0) + r2))*(-omega + omega0 + 1i*r2)^2)"
	),
	sp1=list(	
		func="(aq*((end - off)*pi*cos(end*pi) - exp(aq*(1i*omega - 1i*omega0 + r2))*(end - off)*pi*cos(off*pi) + aq*(1i*omega - 1i*omega0 + r2)*(sin(end*pi) - exp(aq*(1i*omega - 1i*omega0 + r2))*sin(off*pi))))/(exp(aq*(1i*omega - 1i*omega0 + r2))*((end - off)*pi + aq*(omega - omega0 - 1i*r2))*((-end + off)*pi + aq*(omega - omega0 - 1i*r2)))",
		domega0="(aq^2*(-1i*(end - off)*pi*((end - off)^2*pi^2 - aq^2*(-omega + omega0 + 1i*r2)^2 + 2*aq*(1i*(omega - omega0) + r2))*cos(end*pi) - 2*aq*exp(aq*(1i*(omega - omega0) + r2))*(end - off)*pi*(omega - omega0 - 1i*r2)*cos(off*pi) + (1i*(end - off)^2*pi^2 - aq*(-((end - off)^2*pi^2) + aq*(-1i + aq*(omega - omega0 - 1i*r2))*(omega - omega0 - 1i*r2))*(omega - omega0 - 1i*r2))*sin(end*pi) + exp(aq*(1i*(omega - omega0) + r2))*(1i*(end - off)*pi + aq*(omega - omega0 - 1i*r2))*(-1i*aq*omega + 1i*aq*omega0 - end*pi + off*pi - aq*r2)*sin(off*pi)))/(exp(aq*(1i*(omega - omega0) + r2))*((end - off)*pi + aq*(omega - omega0 - 1i*r2))^2*((end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2)",
		dr2="(aq^2*((end - off)*pi*((end - off)^2*pi^2 - aq^2*(-omega + omega0 + 1i*r2)^2 + 2*aq*(1i*(omega - omega0) + r2))*cos(end*pi) - 2i*aq*exp(aq*(1i*(omega - omega0) + r2))*(end - off)*pi*(omega - omega0 - 1i*r2)*cos(off*pi) + (-((end - off)^2*pi^2) + 1i*aq*(omega - omega0 - 1i*r2)*((end - off)^2*pi^2 - aq^2*(-omega + omega0 + 1i*r2)^2 + aq*(1i*(omega - omega0) + r2)))*sin(end*pi) + exp(aq*(1i*(omega - omega0) + r2))*((end - off)^2*pi^2 + aq^2*(-omega + omega0 + 1i*r2)^2)*sin(off*pi)))/(exp(aq*(1i*(omega - omega0) + r2))*((end - off)*pi + aq*(omega - omega0 - 1i*r2))^2*((end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2)"
	),
	sp2=list(
		func="((-1i*aq*(-1 + exp(-(1i*(2*(end - off)*pi + aq*(omega - omega0 - 1i*r2))))))/(exp(2*1i*off*pi)*(2*(end - off)*pi + aq*(omega - omega0 - 1i*r2))) - (1i*aq*exp(2*1i*off*pi)*(-1 + exp(1i*(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2)))))/(2*(-end + off)*pi + aq*(omega - omega0 - 1i*r2)) + (2i*(-1 + exp(-(1i*aq*(omega - omega0 - 1i*r2)))))/(omega - omega0 - 1i*r2))/4",
		domega0="(exp(1i*aq*(-omega + omega0 + 1i*r2))*(aq*((aq*exp(2*1i*end*pi)*(-1i - 2*end*pi + 2*off*pi + aq*(omega - omega0 - 1i*r2)))/(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2 + (aq*(-1i + 2*end*pi - 2*off*pi + aq*(omega - omega0 - 1i*r2)))/(exp(2*1i*end*pi)*(2*(-end + off)*pi + aq*(-omega + omega0 + 1i*r2))^2) - 2/(omega - omega0 - 1i*r2)) + 2i/(-omega + omega0 + 1i*r2)^2) + (1i*aq^2*exp(2*1i*off*pi))/(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2 + (1i*aq^2)/(exp(2*1i*off*pi)*(2*(-end + off)*pi + aq*(-omega + omega0 + 1i*r2))^2) - 2i/(-omega + omega0 + 1i*r2)^2)/4",
		dr2="((2*(1 - exp(-(1i*aq*(omega - omega0 - 1i*r2)))))/(-omega + omega0 + 1i*r2)^2 + (2*aq)/(exp(aq*(1i*(omega - omega0) + r2))*(1i*(omega - omega0) + r2)) + aq^2*(-(exp(2*1i*off*pi)/(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2) - 1/(exp(2*1i*off*pi)*(2*(-end + off)*pi + aq*(-omega + omega0 + 1i*r2))^2) + (exp(1i*(2*end*pi + aq*(-omega + omega0 + 1i*r2)))*(1 - 2i*(end - off)*pi + aq*(1i*(omega - omega0) + r2)))/(2*(end - off)*pi + aq*(-omega + omega0 + 1i*r2))^2 + (exp((-2i)*end*pi - aq*(1i*(omega - omega0) + r2))*(1 + 2i*(end - off)*pi + aq*(1i*(omega - omega0) + r2)))/(2*(-end + off)*pi + aq*(-omega + omega0 + 1i*r2))^2))/4"
	)
)

#lineshapes_parsed <- lapply(lineshapes, function(flist) lapply(flist, function(fvec) parse(text=fvec, keep.source=FALSE)))

#lineshapes_simplified_func <- lapply(lineshapes, function(lineshape_list) compile(findSubexprs(eval(parse(text=paste("quote(", lineshape_list[[1]], ")"))), simplify=TRUE), options=list(optimize=3, suppressUndefined=TRUE)))

#lineshapes_simplified <- lapply(lineshapes, function(lineshape_list) compile(findSubexprs(eval(parse(text=paste("quote(list(", paste(unlist(lineshape_list), collapse=","), "))"))), simplify=TRUE), options=list(optimize=3, suppressUndefined=TRUE)))
#lineshapes_simplified <- lapply(lineshapes, function(lineshape_list) findSubexprs(eval(parse(text=paste("quote(list(", paste(unlist(paste(c("func", "domega0", "dr2"), "=", lineshape_list)), collapse=","), "))"))), simplify=TRUE))

#lineshapes_simplified_func <- lapply(fitnmr:::lineshapes, function(lineshape_list) nlsr::findSubexprs(eval(parse(text=paste("quote(", lineshape_list[[1]], ")"))), simplify=TRUE))
#lineshapes_simplified_func_deriv <- lapply(fitnmr:::lineshapes, function(lineshape_list) nlsr::findSubexprs(eval(parse(text=paste("quote(list(", paste(unlist(paste(c("func", "domega0", "dr2"), "=", lineshape_list)), collapse=","), "))"))), simplify=TRUE))

value_none <- function(omega, omega0, r2, aq, end=NULL, off=NULL) {

	.expr1 <- omega - omega0
	.expr2 <- (0+1i) * r2
	.expr3 <- .expr1 - .expr2
	(0+1i) * (-1 + exp(-((0+1i) * aq * .expr3)))/.expr3
}

value_deriv_none <- function(omega, omega0, r2, aq, end=NULL, off=NULL) {

	.expr1 <- omega - omega0
	.expr2 <- (0+1i) * r2
	.expr3 <- .expr1 - .expr2
	.expr4 <- (0+1i) * .expr1
	.expr5 <- .expr4 + r2
	.expr6 <- aq * .expr5
	.expr7 <- exp(.expr6)
	.expr8 <- -omega
	.expr9 <- .expr8 + omega0
	.expr10 <- .expr9 + .expr2
	.expr11 <- (0+1i) * aq
	.expr12 <- .expr11 * .expr3
	.expr13 <- .expr10^2
	.expr14 <- .expr7 * .expr13
	list(func = (0+1i) * (-1 + exp(-.expr12))/.expr3, domega0 = (0+1i - 
		(0+1i) * .expr7 + aq * .expr10)/.expr14, dr2 = (-1 + 
		.expr7 - .expr12)/.expr14)
}

value_sp1 <- function(omega, omega0, r2, aq, end, off) {

	.expr1 <- end - off
	.expr2 <- (0+1i) * omega
	.expr3 <- (0+1i) * omega0
	.expr4 <- .expr2 - .expr3
	.expr5 <- .expr4 + r2
	.expr6 <- aq * .expr5
	.expr7 <- end * pi
	.expr8 <- exp(.expr6)
	.expr9 <- off * pi
	.expr10 <- .expr1 * pi
	.expr11 <- omega - omega0
	.expr12 <- (0+1i) * r2
	.expr13 <- .expr11 - .expr12
	.expr14 <- aq * .expr13
	aq * (.expr10 * cos(.expr7) - .expr8 * .expr1 * pi * cos(.expr9) + 
		.expr6 * (sin(.expr7) - .expr8 * sin(.expr9)))/(.expr8 * 
		(.expr10 + .expr14) * ((-end + off) * pi + .expr14))
}

value_deriv_sp1 <- function(omega, omega0, r2, aq, end, off) {

	.expr1 <- end - off
	.expr2 <- (0+1i) * omega
	.expr3 <- (0+1i) * omega0
	.expr4 <- .expr2 - .expr3
	.expr5 <- .expr4 + r2
	.expr6 <- aq * .expr5
	.expr7 <- end * pi
	.expr8 <- exp(.expr6)
	.expr9 <- off * pi
	.expr10 <- .expr1 * pi
	.expr11 <- omega - omega0
	.expr12 <- (0+1i) * r2
	.expr13 <- .expr11 - .expr12
	.expr14 <- aq * .expr13
	.expr15 <- aq^2
	.expr16 <- cos(.expr7)
	.expr17 <- 2 * aq
	.expr18 <- (0+1i) * .expr11
	.expr19 <- .expr18 + r2
	.expr20 <- cos(.expr9)
	.expr21 <- .expr1^2
	.expr22 <- pi^2
	.expr23 <- .expr21 * .expr22
	.expr24 <- -(0+1i)
	.expr25 <- sin(.expr7)
	.expr26 <- aq * .expr19
	.expr27 <- exp(.expr26)
	.expr28 <- sin(.expr9)
	.expr29 <- .expr10 + .expr14
	.expr30 <- -omega
	.expr31 <- .expr30 + omega0
	.expr32 <- .expr31 + .expr12
	.expr33 <- .expr32^2
	.expr34 <- .expr15 * .expr33
	.expr35 <- .expr23 - .expr34
	.expr36 <- .expr17 * .expr19
	.expr37 <- .expr35 + .expr36
	.expr38 <- -.expr23
	.expr39 <- (0+1i) * aq
	.expr40 <- .expr29^2
	.expr41 <- .expr27 * .expr40
	.expr42 <- aq * .expr32
	.expr43 <- .expr10 + .expr42
	.expr44 <- .expr43^2
	.expr45 <- .expr41 * .expr44
	list(func = aq * (.expr10 * .expr16 - .expr8 * .expr1 * pi * 
		.expr20 + .expr6 * (.expr25 - .expr8 * .expr28))/(.expr8 * 
		.expr29 * ((-end + off) * pi + .expr14)), domega0 = .expr15 * 
		(.expr24 * .expr1 * pi * .expr37 * .expr16 - .expr17 * 
			.expr27 * .expr1 * pi * .expr13 * .expr20 + ((0+1i) * 
			.expr21 * .expr22 - aq * (.expr38 + aq * (.expr24 + 
			.expr14) * .expr13) * .expr13) * .expr25 + .expr27 * 
			((0+1i) * .expr1 * pi + .expr14) * (.expr24 * aq * 
			omega + .expr39 * omega0 - .expr7 + .expr9 - aq * 
			r2) * .expr28)/.expr45, dr2 = .expr15 * (.expr10 * 
		.expr37 * .expr16 - (0+2i) * aq * .expr27 * .expr1 * 
		pi * .expr13 * .expr20 + (.expr38 + .expr39 * .expr13 * 
		(.expr35 + .expr26)) * .expr25 + .expr27 * (.expr23 + 
		.expr34) * .expr28)/.expr45)
}

value_sp2 <- function(omega, omega0, r2, aq, end, off) {

	.expr1 <- end - off
	.expr2 <- 2 * .expr1
	.expr3 <- .expr2 * pi
	.expr4 <- omega - omega0
	.expr5 <- (0+1i) * r2
	.expr6 <- .expr4 - .expr5
	.expr7 <- aq * .expr6
	.expr8 <- .expr3 + .expr7
	.expr9 <- 2 * (0+1i)
	.expr10 <- .expr9 * off
	.expr11 <- .expr10 * pi
	.expr12 <- exp(.expr11)
	.expr13 <- (0+1i) * aq
	(-(0+1i) * aq * (-1 + exp(-((0+1i) * .expr8)))/(.expr12 * 
		.expr8) - .expr13 * .expr12 * (-1 + exp((0+1i) * (.expr3 + 
		aq * (-omega + omega0 + .expr5))))/(2 * (-end + off) * 
		pi + .expr7) + (0+2i) * (-1 + exp(-(.expr13 * .expr6)))/.expr6)/4
}

value_deriv_sp2 <- function(omega, omega0, r2, aq, end, off) {

	.expr1 <- end - off
	.expr2 <- 2 * .expr1
	.expr3 <- .expr2 * pi
	.expr4 <- omega - omega0
	.expr5 <- (0+1i) * r2
	.expr6 <- .expr4 - .expr5
	.expr7 <- aq * .expr6
	.expr8 <- .expr3 + .expr7
	.expr9 <- 2 * (0+1i)
	.expr10 <- .expr9 * off
	.expr11 <- .expr10 * pi
	.expr12 <- exp(.expr11)
	.expr13 <- (0+1i) * aq
	.expr14 <- -omega
	.expr15 <- .expr14 + omega0
	.expr16 <- .expr15 + .expr5
	.expr17 <- -(0+1i)
	.expr18 <- aq * .expr16
	.expr19 <- .expr3 + .expr18
	.expr20 <- 2 * end
	.expr21 <- .expr20 * pi
	.expr22 <- 2 * off
	.expr23 <- .expr22 * pi
	.expr24 <- .expr9 * end
	.expr25 <- .expr24 * pi
	.expr26 <- exp(.expr25)
	.expr27 <- -end
	.expr28 <- .expr27 + off
	.expr29 <- 2 * .expr28
	.expr30 <- .expr29 * pi
	.expr31 <- .expr19^2
	.expr32 <- aq^2
	.expr33 <- (0+1i) * .expr32
	.expr34 <- .expr30 + .expr18
	.expr35 <- .expr34^2
	.expr36 <- .expr16^2
	.expr37 <- (0+2i)/.expr36
	.expr38 <- .expr13 * .expr6
	.expr39 <- -.expr38
	.expr40 <- exp(.expr39)
	.expr41 <- (0+1i) * .expr4
	.expr42 <- .expr41 + r2
	.expr43 <- .expr12 * .expr35
	.expr44 <- aq * .expr42
	.expr45 <- (0+2i) * .expr1
	.expr46 <- .expr45 * pi
	list(func = (.expr17 * aq * (-1 + exp(-((0+1i) * .expr8)))/(.expr12 * 
		.expr8) - .expr13 * .expr12 * (-1 + exp((0+1i) * .expr19))/(.expr30 + 
		.expr7) + (0+2i) * (-1 + .expr40)/.expr6)/4, domega0 = (exp(.expr13 * 
		.expr16) * (aq * (aq * .expr26 * (.expr17 - .expr21 + 
		.expr23 + .expr7)/.expr31 + aq * (.expr17 + .expr21 - 
		.expr23 + .expr7)/(.expr26 * .expr35) - 2/.expr6) + .expr37) + 
		.expr33 * .expr12/.expr31 + .expr33/.expr43 - .expr37)/4, 
		dr2 = (2 * (1 - .expr40)/.expr36 + 2 * aq/(exp(.expr44) * 
			.expr42) + .expr32 * (-(.expr12/.expr31) - 1/.expr43 + 
			exp((0+1i) * (.expr21 + .expr18)) * (1 - .expr46 + 
				.expr44)/.expr31 + exp(-(0+2i) * end * pi - .expr44) * 
			(1 + .expr46 + .expr44)/.expr35))/4)
}