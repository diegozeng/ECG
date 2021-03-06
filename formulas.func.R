formulas <- function(){
  f1 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24+V25+V26+V27+V28+V29+V30+V31+V32+V33+V34+V35+V36+V37+V38+V39+V40+V41+V42+V43+V44
  f2 <- V45 ~ V4+V5+V6+V7+V8+V9+V10#RR
  f3 <- V45 ~ V11+V12+V13+V14+V15+V16+V17#QT
  f4 <- V45 ~ V18+V19+V20+V21+V22+V23+V24#QT/RR
  f5 <- V45 ~ V25+V26+V27+V28+V29+V30+V31+V32+V33+V34+V35+V36+V37+V38+V39+V40+V41+V42+V43+V44#Reg
  f6 <- V45 ~ V25+V26#Reg_lin
  f7 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24#RR+QT+QT/RR
  f8 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17#RR+QT
  f9 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V18+V19+V20+V21+V22+V23+V24#RR+QT/RR
  f10 <- V45 ~ V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24#QT+QT/RR
  f11 <- V45 ~ V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24+V25+V26#RR+QT+QT/RR+Reg_Lin
  
  c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11)
}