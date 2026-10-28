# sumChains
#' Calculate total number of carbons and double bounds of lipid chains
#'
#' Given the structure of a lipid specie, it sums up the chains.
#'
#' @param chains character value with the configuration of the chains separated
#' by a white space
#' @param n number of chains
#'
#' @return Character value indicating the total number of carbons and double
#' bounds
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
sumChains <- function(chains, n) {
  cb <- unlist(strsplit(chains, "[: ]"))
  if (n == 1) {
    sum_c <- as.numeric(cb[1])
    sum_b <- as.numeric(cb[2])
  }
  if (n == 2) {
    sum_c <- as.numeric(cb[1]) + as.numeric(cb[3])
    sum_b <- as.numeric(cb[2]) + as.numeric(cb[4])
  } else if (n == 3) {
    sum_c <- as.numeric(cb[1]) + as.numeric(cb[3]) + as.numeric(cb[5])
    sum_b <- as.numeric(cb[2]) + as.numeric(cb[4]) + as.numeric(cb[6])
  } else if (n == 4) {
    sum_c <- as.numeric(cb[1]) + as.numeric(cb[3]) + as.numeric(cb[5]) + as.numeric(cb[7])
    sum_b <- as.numeric(cb[2]) + as.numeric(cb[4]) + as.numeric(cb[6]) + as.numeric(cb[8])
  }
  total <- paste(c(sum_c, sum_b), collapse=":")
  return(total)
}

# MassSph_d
#' Calculate formula and mass of sphingoid bases
#'
#' Calculate formula and mass of sphingoid bases
#'
#' @param Sph character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSph_d <- function(Sph) {
  cb <- unlist(strsplit(Sph, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2 + 3 - as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO2"),
                   collapse="")
  mass <- C*12 + H*1.007825 + 14.00305 + 2*15.9949
  return(c(formula, mass))
}

# MassSph_t
#' Calculate formula and mass of sphingoid bases
#'
#' Calculate formula and mass of sphingoid bases
#'
#' @param Sph character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSph_t <- function(Sph) {
  cb <- unlist(strsplit(Sph, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2 + 5 - as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO3"),
                   collapse="")
  mass <- C*12 + H*1.007825 + 14.00305 + 3*15.9949
  return(c(formula, mass))
}

# MassSphP_d
#' Calculate formula and mass of sphingoid phosphate bases
#'
#' Calculate formula and mass of sphingoid phosphate bases
#'
#' @param SphP character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSphP_d <- function(SphP) {
  cb <- unlist(strsplit(SphP, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2 + 4 - as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO5P"),
                   collapse="")
  mass <- C*12 + H*1.007825 + 14.00305 + 5*15.9949 + 30.97393
  return(c(formula, mass))
}

# MassSphP_t
#' Calculate formula and mass of sphingoid phosphate bases
#'
#' Calculate formula and mass of sphingoid phosphate bases
#'
#' @param SphP character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSphP_t <- function(SphP) {
  cb <- unlist(strsplit(SphP, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2 + 6 - as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12 + H*1.007825 + 14.00305 + 6*15.9949 + 30.97393
  return(c(formula, mass))
}

# MassCer_d
#' Calculate formula and mass of Cer_d
#'
#' Calculate formula and mass of Cer_d
#'
#' @param cer character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCer_d <- function(cer) {
  cb <- unlist(strsplit(cer, "[: ]"))
  C <- as.numeric(cb[1])
  H <- C*2 + 1 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO3")
  mass <- C*12 + H*1.007825 + 14.00305 + 3*15.9949
  return(c(formula, mass))
}

# MassCer_t
#' Calculate formula and mass of Cer_t
#'
#' Calculate formula and mass of Cer_t
#'
#' @param cer character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCer_t <- function(cer) {
  cb <- unlist(strsplit(cer, "[: ]"))
  C <- as.numeric(cb[1])
  H <- C*2 + 3 - as.numeric(cb[2])*2   # +2H vs d
  formula <- paste0("C", C, "H", H, "NO4")  # +O
  mass <- C*12 + H*1.007825 + 14.00305 + 4*15.9949
  return(c(formula, mass))
}

# MassHexCer_d
#' Calculate formula and mass of HexCer_d
#'
#' Calculate formula and mass of HexCer_d
#'
#' @param hexcer character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassHexCer_d <- function(hexcer) {
  cb <- unlist(strsplit(hexcer, "[: ]"))
  C <- as.numeric(cb[1]) + 6
  H <- as.numeric(cb[1])*2 + 11 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO8")
  mass <- C*12 + H*1.007825 + 14.00305 + 8*15.9949
  return(c(formula, mass))
}

# MassHexCer_t
#' Calculate formula and mass of HexCer_t
#'
#' Calculate formula and mass of HexCer_t
#'
#' @param hexcer character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassHexCer_t <- function(hexcer) {
  cb <- unlist(strsplit(hexcer, "[: ]"))
  C <- as.numeric(cb[1]) + 6
  H <- as.numeric(cb[1])*2 + 13 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO9")
  mass <- C*12 + H*1.007825 + 14.00305 + 9*15.9949
  return(c(formula, mass))
}

# MassAHexCer_d
#' Calculate formula and mass of AHexCer_d
#'
#' Calculate formula and mass of AHexCer_d
#'
#' @param ahexcer character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassAHexCer_d <- function(ahexcer) {
  cb <- unlist(strsplit(ahexcer, "[: ]"))
  C <- as.numeric(cb[1]) + 6
  H <- as.numeric(cb[1])*2 + 9 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO9")
  mass <- C*12 + H*1.007825 + 14.00305 + 9*15.9949
  return(c(formula, mass))
}

# MassAHexCer_t
#' Calculate formula and mass of AHexCer_t
#'
#' Calculate formula and mass of AHexCer_t
#'
#' @param ahexcer character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassAHexCer_t <- function(ahexcer) {
  cb <- unlist(strsplit(ahexcer, "[: ]"))
  C <- as.numeric(cb[1]) + 6
  H <- as.numeric(cb[1])*2 + 11 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO10")
  mass <- C*12 + H*1.007825 + 14.00305 + 10*15.9949
  return(c(formula, mass))
}

# MassSHexCer_d
#' Calculate formula and mass of SHexCer_d
#'
#' Calculate formula and mass of SHexCer_d
#'
#' @param shexcer character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSHexCer_d <- function(shexcer) {
  cb <- unlist(strsplit(shexcer, "[: ]"))
  C <- as.numeric(cb[1]) + 6
  H <- as.numeric(cb[1])*2 + 11 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO11S")
  mass <- C*12 + H*1.007825 + 14.00305 + 11*15.9949 + 31.97207
  return(c(formula, mass))
}

# MassSHexCer_t
#' Calculate formula and mass of SHexCer_t
#'
#' Calculate formula and mass of SHexCer_t
#'
#' @param shexcer character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSHexCer_t <- function(shexcer) {
  cb <- unlist(strsplit(shexcer, "[: ]"))
  C <- as.numeric(cb[1]) + 6
  H <- as.numeric(cb[1])*2 + 13 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO12S")
  mass <- C*12 + H*1.007825 + 14.00305 + 12*15.9949 + 31.97207
  return(c(formula, mass))
}

# MassCerP_d
#' Calculate formula and mass of ceramides phosphate
#'
#' Calculate formula and mass of ceramides phosphate
#'
#' @param cerP character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCerP_d <- function(cerP) {
  cb <- unlist(strsplit(cerP, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2 + 2 - as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12 + H*1.007825 + 14.00305 + 6*15.9949 + 30.97393
  return(c(formula, mass))
}

# MassCerP_t
#' Calculate formula and mass of ceramides phosphate
#'
#' Calculate formula and mass of ceramides phosphate
#'
#' @param cerP character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCerP_t <- function(cerP) {
  cb <- unlist(strsplit(cerP, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2 + 4 - as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12 + H*1.007825 + 14.00305 + 7*15.9949 + 30.97393
  return(c(formula, mass))
}

# MassAcylCer_d
#' Calculate formula and mass of acylceramides
#'
#' Calculate formula and mass of acylceramides
#'
#' @param acylcer character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassAcylCer_d <- function(acylcer) {
  cb <- unlist(strsplit(acylcer, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2 - as.numeric(cb[2])*2 - 1
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO4"),
                   collapse="")
  mass <- C*12 + H*1.007825 + 14.00305 + 4*15.9949
  return(c(formula, mass))
}

# MassAcylCer_t
#' Calculate formula and mass of acylceramides
#'
#' Calculate formula and mass of acylceramides
#'
#' @param acylcer character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassAcylCer_t <- function(acylcer) {
  cb <- unlist(strsplit(acylcer, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2 + 2 - as.numeric(cb[2])*2 - 1
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO5"),
                   collapse="")
  mass <- C*12 + H*1.007825 + 14.00305 + 5*15.9949
  return(c(formula, mass))
}

# MassSM_d
#' Calculate formula and mass of SM_d
#'
#' Calculate formula and mass of SM_d
#'
#' @param SM character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSM_d <- function(SM) {
  cb <- unlist(strsplit(SM, "[: ]"))
  C <- as.numeric(cb[1]) + 5
  H <- as.numeric(cb[1]) * 2 + 13 - as.numeric(cb[2]) * 2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "N2O6P"),
                   collapse = "")
  mass <- C*12 + H*1.007825 + 2*14.00305 + 6*15.9949 + 30.97393
  return(c(formula, mass))
}

# MassSM_t
#' Calculate formula and mass of SM_t
#'
#' Calculate formula and mass of SM_t
#'
#' @param SM character value indicating total number of carbons and double
#' bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassSM_t <- function(SM) {
  cb <- unlist(strsplit(SM, "[: ]"))
  C <- as.numeric(cb[1]) + 5
  H <- as.numeric(cb[1]) * 2 + 15 - as.numeric(cb[2]) * 2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "N2O7P"),
                   collapse = "")  # +1 O
  mass <- C*12 + H*1.007825 + 2*14.00305 + 7*15.9949 + 30.97393
  return(c(formula, mass))
}

# MassCarnitines
#' Calculate formula and mass of carnitines
#'
#' Calculate formula and mass of carnitines
#'
#' @param carnitine character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCarnitines <- function(carnitine) {
  cb <- unlist(strsplit(carnitine, "[: ]"))
  C <- as.numeric(cb[1])+7
  H <- as.numeric(cb[1])*2+13-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO4"),
                   collapse="")
  mass <- C*12+H*1.007825++14.00305+4*15.9949
  return(c(formula, mass))
}

# MassCE
#' Calculate formula and mass of cholesterol esthers
#'
#' Calculate formula and mass of cholesterol esthers
#'
#' @param CE character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCE <- function(CE) {
  cb <- unlist(strsplit(CE, "[: ]"))
  C <- as.numeric(cb[1])+27
  H <- as.numeric(cb[1])*2+44-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O2"),
                   collapse="")
  mass <- C*12+H*1.007825+2*15.9949
  return(c(formula, mass))
}

# MassFA
#' Calculate formula and mass of fatty acids
#'
#' Calculate formula and mass of fatty acids
#'
#' @param FA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassFA <- function(FA) {
  cb <- unlist(strsplit(FA, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O2"),
                   collapse="")
  mass <- C*12+H*1.007825+2*15.9949
  return(c(formula, mass))
}

# MassHFA
#' Calculate formula and mass of hydroxi fatty acids
#'
#' Calculate formula and mass of hydroxi fatty acids
#'
#' @param HFA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassHFA <- function(HFA) {
  cb <- unlist(strsplit(HFA, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O3"),
                   collapse="")
  mass <- C*12+H*1.007825+3*15.9949
  return(c(formula, mass))
}

# MassFAHFA
#' Calculate formula and mass of FAHFA
#'
#' Calculate formula and mass of FAHFA
#'
#' @param FAHFA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassFAHFA <- function(FAHFA) {
  cb <- unlist(strsplit(FAHFA, "[: ]"))
  C <- as.numeric(cb[1])
  H <- as.numeric(cb[1])*2-2-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O4"),
                   collapse="")
  mass <- C*12+H*1.007825+4*15.9949
  return(c(formula, mass))
}

# MassLysoPA
#' Calculate formula and mass of LPA
#'
#' Calculate formula and mass of LPA
#'
#' @param LPA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPA <- function(LPA) {
  cb <- unlist(strsplit(LPA, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+7-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H),
                     "O7P"), collapse="")
  mass <- C*12+H*1.007825+7*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPAo
#' Calculate formula and mass of LPAo
#'
#' Calculate formula and mass of LPAo
#'
#' @param LPAo character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPAo <- function(LPAo) {
  cb <- unlist(strsplit(LPAo, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+9-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H),
                     "O6P"), collapse="")
  mass <- C*12+H*1.007825+6*15.9949+30.97393
  return(c(formula, mass))
}

# MassPA
#' Calculate formula and mass of PA
#'
#' Calculate formula and mass of PA
#'
#' @param PA character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPA <- function(PA) {
  cb <- unlist(strsplit(PA, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+5-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O8P"),
                   collapse="")
  mass <- C*12+H*1.007825+8*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPE
#' Calculate formula and mass of LPE
#'
#' Calculate formula and mass of LPE
#'
#' @param LPE character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPE <- function(LPE) {
  cb <- unlist(strsplit(LPE, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+12-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+7*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPEo
#' Calculate formula and mass of LPEo
#'
#' Calculate formula and mass of LPEo
#'
#' @param LPEo character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPEo <- function(LPEo) {
  cb <- unlist(strsplit(LPEo, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+14-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+6*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPEp
#' Calculate formula and mass of LPEp
#'
#' Calculate formula and mass of LPEp
#'
#' @param LPEp character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPEp <- function(LPEp) {
  cb <- unlist(strsplit(LPEp, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+12-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+6*15.9949+30.97393
  return(c(formula, mass))
}

# MassPE
#' Calculate formula and mass of PE
#'
#' Calculate formula and mass of PE
#'
#' @param PE character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPE <- function(PE) {
  cb <- unlist(strsplit(PE, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+10-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO8P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+8*15.9949+30.97393
  return(c(formula, mass))
}

# MassPEo
#' Calculate formula and mass of plasmanyl PE
#'
#' Calculate formula and mass of plasmanyl PE
#'
#' @param PEo character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPEo <- function(PEo) {
  cb <- unlist(strsplit(PEo, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+12-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+7*15.9949+30.97393
  return(c(formula, mass))
}

# MassPEp
#' Calculate formula and mass of plasmenyl PE
#'
#' Calculate formula and mass of plasmenyl PE
#'
#' @param PEp character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPEp <- function(PEp) {
  cb <- unlist(strsplit(PEp, "[: ]"))
  C <- as.numeric(cb[1])+5
  H <- as.numeric(cb[1])*2+10-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+14.00305+7*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPG
#' Calculate formula and mass of LPG
#'
#' Calculate formula and mass of LPG
#'
#' @param LPG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPG <- function(LPG) {
  cb <- unlist(strsplit(LPG, "[: ]"))
  C <- as.numeric(cb[1])+6
  H <- as.numeric(cb[1])*2+13-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H),
                     "O9P"), collapse="")
  mass <- C*12+H*1.007825+9*15.9949+30.97393
  return(c(formula, mass))
}

# MassPG
#' Calculate formula and mass of PG
#'
#' Calculate formula and mass of PG
#'
#' @param PG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPG <- function(PG) {
  cb <- unlist(strsplit(PG, "[: ]"))
  C <- as.numeric(cb[1])+6
  H <- as.numeric(cb[1])*2+11-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H),
                     "O10P"), collapse="")
  mass <- C*12+H*1.007825+10*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPI
#' Calculate formula and mass of LPI
#'
#' Calculate formula and mass of LPI
#'
#' @param LPI character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPI <- function(LPI) {
  cb <- unlist(strsplit(LPI, "[: ]"))
  C <- as.numeric(cb[1])+9
  H <- as.numeric(cb[1])*2+17-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O12P"),
                   collapse="")
  mass <- C*12+H*1.007825+12*15.9949+30.97393
  return(c(formula, mass))
}

# MassPI
#' Calculate formula and mass of PI
#'
#' Calculate formula and mass of PI
#'
#' @param PI character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPI <- function(PI) {
  cb <- unlist(strsplit(PI, "[: ]"))
  C <- as.numeric(cb[1])+9
  H <- as.numeric(cb[1])*2+15-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O13P"),
                   collapse="")
  mass <- C*12+H*1.007825+13*15.9949+30.97393
  return(c(formula, mass))
}

# MassLysoPS
#' Calculate formula and mass of LysoPS
#'
#' Calculate formula and mass of LysoPS
#'
#' @param LPS character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPS <- function(LPS) {
  cb <- unlist(strsplit(LPS, "[: ]"))
  C <- as.numeric(cb[1])+6
  H <- as.numeric(cb[1])*2+12-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO9P"),
                   collapse="")
  mass <- C*12+H*1.007825+9*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassPS
#' Calculate formula and mass of PS
#'
#' Calculate formula and mass of PS
#'
#' @param PS character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPS <- function(PS) {
  cb <- unlist(strsplit(PS, "[: ]"))
  C <- as.numeric(cb[1])+6
  H <- as.numeric(cb[1])*2+10-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO10P"),
                   collapse="")
  mass <- C*12+H*1.007825+10*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassLysoPC
#' Calculate formula and mass of LysoPC
#'
#' Calculate formula and mass of LysoPC
#'
#' @param LPC character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPC <- function(LPC) {
  cb <- unlist(strsplit(LPC, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+18-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+7*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassLysoPCp
#' Calculate formula and mass of LysoPCp
#'
#' Calculate formula and mass of LysoPCp
#'
#' @param LPCp character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPCp <- function(LPCp) {
  cb <- unlist(strsplit(LPCp, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+18-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12+H*1.007825+6*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassLysoPCo
#' Calculate formula and mass of LysoPCo
#'
#' Calculate formula and mass of LysoPCo
#'
#' @param LPCo character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLysoPCo <- function(LPCo) {
  cb <- unlist(strsplit(LPCo, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+20-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO6P"),
                   collapse="")
  mass <- C*12+H*1.007825+6*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassPC
#' Calculate formula and mass of PC
#'
#' Calculate formula and mass of PC
#'
#' @param PC character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPC <- function(PC) {
  cb <- unlist(strsplit(PC, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+16-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO8P"),
                   collapse="")
  mass <- C*12+H*1.007825+8*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassPCp
#' Calculate formula and mass of PCp
#'
#' Calculate formula and mass of PCp
#'
#' @param PCp character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPCp <- function(PCp) {
  cb <- unlist(strsplit(PCp, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+16-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+7*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassPCo
#' Calculate formula and mass of PCo
#'
#' Calculate formula and mass of PCo
#'
#' @param PCo character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassPCo <- function(PCo) {
  cb <- unlist(strsplit(PCo, "[: ]"))
  C <- as.numeric(cb[1])+8
  H <- as.numeric(cb[1])*2+18-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "NO7P"),
                   collapse="")
  mass <- C*12+H*1.007825+7*15.9949+30.97393+14.00305
  return(c(formula, mass))
}

# MassMG
#' Calculate formula and mass of MG
#'
#' Calculate formula and mass of MG
#'
#' @param MG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassMG <- function(MG) {
  cb <- unlist(strsplit(MG, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+6-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O4"),
                   collapse="")
  mass <- C*12+H*1.007825+4*15.9949
  return(c(formula, mass))
}

# MassDG
#' Calculate formula and mass of DG
#'
#' Calculate formula and mass of DG
#'
#' @param DG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassDG <- function(DG) {
  cb <- unlist(strsplit(DG, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+4-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O5"),
                   collapse="")
  mass <- C*12+H*1.007825+5*15.9949
  return(c(formula, mass))
}

# MassTG
#' Calculate formula and mass of TG
#'
#' Calculate formula and mass of TG
#'
#' @param TG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassTG <- function(TG) {
  cb <- unlist(strsplit(TG, "[: ]"))
  C <- as.numeric(cb[1])+3
  H <- as.numeric(cb[1])*2+2-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O6"),
                   collapse="")
  mass <- C*12+H*1.007825+6*15.9949
  return(c(formula, mass))
}

# MassMGDG
#' Calculate formula and mass of MGDG
#'
#' Calculate formula and mass of MGDG
#'
#' @param MGDG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassMGDG <- function(MGDG) {
  cb <- unlist(strsplit(MGDG, "[: ]"))
  C <- as.numeric(cb[1]) + 6
  H <- as.numeric(cb[1])*2 + 10 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "O10")
  mass <- C*12 + H*1.007825 + 10*15.9949
  return(c(formula, mass))
}

# MassDGDG
#' Calculate formula and mass of DGDG
#'
#' Calculate formula and mass of DGDG
#'
#' @param DGDG character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassDGDG <- function(DGDG) {
  cb <- unlist(strsplit(DGDG, "[: ]"))
  C <- as.numeric(cb[1]) + 12
  H <- as.numeric(cb[1])*2 + 20 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "O15")
  mass <- C*12 + H*1.007825 + 15*15.9949
  return(c(formula, mass))
}

# MassBMP
#' Calculate formula and mass of BMP
#'
#' Calculate formula and mass of BMP
#'
#' @param BMP character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassBMP <- function(BMP) {
  cb <- unlist(strsplit(BMP, "[: ]"))
  C <- as.numeric(cb[1]) + 6
  H <- as.numeric(cb[1])*2 + 8 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "O10P")
  mass <- C*12 + H*1.007825 + 10*15.9949 + 30.97393
  return(c(formula, mass))
}

# MassCL
#' Calculate formula and mass of CL
#'
#' Calculate formula and mass of CL
#'
#' @param CL character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassCL <- function(CL) {
  cb <- unlist(strsplit(CL, "[: ]"))
  C <- as.numeric(cb[1])+9
  H <- as.numeric(cb[1])*2+14-as.numeric(cb[2])*2
  formula <- paste(c("C", as.character(C), "H", as.character(H), "O17P2"),
                   collapse="")
  mass <- C*12+H*1.007825+17*15.9949+2*30.97393
  return(c(formula, mass))
}

# MassNAGly
#' Calculate formula and mass of NAGly
#'
#' Calculate formula and mass of NAGly
#'
#' @param NAGly character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassNAGly <- function(NAGly) {
  cb <- unlist(strsplit(NAGly, "[: ]"))
  C <- as.numeric(cb[1]) + 2
  H <- as.numeric(cb[1])*2 + 3 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO3")
  mass <- C*12 + H*1.007825 + 14.00305 + 3*15.9949
  return(c(formula, mass))
}

# MassNAGlySer
#' Calculate formula and mass of NAGlySer
#'
#' Calculate formula and mass of NAGlySer
#'
#' @param NAGlySer character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassNAGlySer <- function(NAGlySer) {
  cb <- unlist(strsplit(NAGlySer, "[: ]"))
  C <- as.numeric(cb[1]) + 3
  H <- as.numeric(cb[1])*2 + 5 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO4")
  mass <- C*12 + H*1.007825 + 14.00305 + 4*15.9949
  return(c(formula, mass))
}

# MassNAOrn
#' Calculate formula and mass of NAOrn
#'
#' Calculate formula and mass of NAOrn
#'
#' @param NAOrn character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassNAOrn <- function(NAOrn) {
  cb <- unlist(strsplit(NAOrn, "[: ]"))
  C <- as.numeric(cb[1]) + 5
  H <- as.numeric(cb[1])*2 + 10 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "N2O2")
  mass <- C*12 + H*1.007825 + 2*14.00305 + 2*15.9949
  return(c(formula, mass))
}

# MassNAE
#' Calculate formula and mass of NAE
#'
#' Calculate formula and mass of NAE
#'
#' @param NAE character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassNAE <- function(NAE) {
  cb <- unlist(strsplit(NAE, "[: ]"))
  C <- as.numeric(cb[1]) + 2
  H <- as.numeric(cb[1])*2 + 5 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO2")
  mass <- C*12 + H*1.007825 + 14.00305 + 2*15.9949
  return(c(formula, mass))
}

# MassLNAPE
#' Calculate formula and mass of LNAPE
#'
#' Calculate formula and mass of LNAPE
#'
#' @param LNAPE character value indicating total number of carbons and
#' double bounds
#'
#' @return vector containing formula and mass
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
MassLNAPE <- function(LNAPE) {
  cb <- unlist(strsplit(LNAPE, "[: ]"))
  C <- as.numeric(cb[1]) + 7
  H <- as.numeric(cb[1])*2 + 13 - as.numeric(cb[2])*2
  formula <- paste0("C", C, "H", H, "NO8P")
  mass <- C*12 + H*1.007825 + 14.00305 + 8*15.9949 + 30.97393
  return(c(formula, mass))
}

# dbSphingolipids
#' Creation of a database for Cer, CerP, GlcCer and SM
#'
#' Creation of a database for Cer, CerP, GlcCer and SM
#'
#' @param chains character vector indicating the FAs to be employed
#' @param chains2 character vector indicating the sphingoid bases to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbSphingolipids <- function(chains, chains2, lipid) {
  comb <- vector()
  if (lipid %in% c("AcylCer_d", "AcylCer_t")){
    comb1 <- unique(utils::combn(chains, 2, simplify=F))
    chains <- unique(unlist(lapply(comb1, paste, collapse=" ")))
    for (i in chains2){
      comb <- append(comb, sapply(chains, paste, i, collapse=" "))
    }
    total  <- unique(unlist(lapply(comb, sumChains, n = 3)))
  } else {
    for (i in chains2){
      comb <- append(comb, sapply(chains, paste, i, collapse=" "))
    }
    comb <- unique(comb)
    total  <- unique(unlist(lapply(comb, sumChains, n = 2)))
  }
  if (lipid == "Cer_d"){
    fm <- lapply(total, MassCer_d)
  } else if (lipid == "Cer_t"){
    fm <- lapply(total, MassCer_t)
  } else if (lipid == "HexCer_d"){
    fm <- lapply(total, MassHexCer_d)
  } else if (lipid == "HexCer_t"){
    fm <- lapply(total, MassHexCer_t)
  } else if (lipid == "AHexCer_d"){
    fm <- lapply(total, MassAHexCer_d)
  } else if (lipid == "AHexCer_t"){
    fm <- lapply(total, MassAHexCer_t)
  } else if (lipid == "SHexCer_d"){
    fm <- lapply(total, MassSHexCer_d)
  } else if (lipid == "SHexCer_t"){
    fm <- lapply(total, MassSHexCer_t)
  } else if (lipid == "CerP_d"){
    fm <- lapply(total, MassCerP_d)
  } else if (lipid == "CerP_t"){
    fm <- lapply(total, MassCerP_t)
  } else if (lipid == "SM_d"){
    fm <- lapply(total, MassSM_d)
  } else if (lipid == "SM_t"){
    fm <- lapply(total, MassSM_t)
  } else if (lipid == "AcylCer_d"){
    fm <- lapply(total, MassAcylCer_d)
  } else if (lipid == "AcylCer_t"){
    fm <- lapply(total, MassAcylCer_t)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=total,
                   Mass=as.numeric(unlist(lapply(fm, "[[", 2))), stringsAsFactors = F)
  return(db)
}

# dbOneChain
#' Creation of a database for Carnitines, CE, FA, HFA, LPL, MG, sphingoid
#' bases and sphingoid bases phosphate.
#'
#' Creation of a database for Carnitines, CE, FA, HFA, LPL, MG, sphingoid
#' bases and sphingoid bases phosphate.
#'
#' @param chains character vector indicating the FAs to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbOneChain <- function(chains, lipid) {
  fas <- unique(chains)
  if (lipid == "Carnitine") {
    fm <- lapply(fas, MassCarnitines)
  } else if (lipid == "CE") {
    fm <- lapply(fas, MassCE)
  } else if (lipid == "FA") {
    fm <- lapply(fas, MassFA)
  } else if (lipid == "HFA") {
    fm <- lapply(fas, MassHFA)
  } else if(lipid == "LPA") {
    fm <- lapply(fas, MassLysoPA)
  } else if(lipid == "LPE") {
    fm <- lapply(fas, MassLysoPE)
  } else if(lipid == "LPG") {
    fm <- lapply(fas, MassLysoPG)
  } else if(lipid == "LPI") {
    fm <- lapply(fas, MassLysoPI)
  } else if(lipid == "LPS") {
    fm <- lapply(fas, MassLysoPS)
  } else if(lipid == "LPC") {
    fm <- lapply(fas, MassLysoPC)
  } else if(lipid == "MG") {
    fm <- lapply(fas, MassMG)
  } else if(lipid == "Sph_d") {
    fm <- lapply(fas, MassSph_d)
  } else if(lipid == "Sph_t") {
    fm <- lapply(fas, MassSph_t)
  } else if(lipid == "SphP_d") {
    fm <- lapply(fas, MassSphP_d)
  } else if(lipid == "SphP_t") {
    fm <- lapply(fas, MassSphP_t)
  } else if(lipid == "LPC-O") {
    fm <- lapply(fas, MassLysoPCo)
  } else if(lipid == "LPC-P") {
    fm <- lapply(fas, MassLysoPCp)
  } else if(lipid == "LPE-O") {
    fm <- lapply(fas, MassLysoPEo)
  } else if(lipid == "LPE-P") {
    fm <- lapply(fas, MassLysoPEp)
  } else if(lipid == "LPA-O") {
    fm <- lapply(fas, MassLysoPAo)
  } else if(lipid == "NAGly") {
    fm <- lapply(fas, MassNAGly)
  } else if(lipid == "NAGlySer") {
    fm <- lapply(fas, MassNAGlySer)
  } else if(lipid == "NAOrn") {
    fm <- lapply(fas, MassNAOrn)
  } else if(lipid == "NAE") {
    fm <- lapply(fas, MassNAE)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=fas,
                   Mass=as.numeric(as.vector(unlist(lapply(fm, "[[", 2)))),
                   stringsAsFactors = F)
  return(db)
}

# dbTwoChains
#' Creation of a database for FAHFA, DG and PL.
#'
#' Creation of a database for FAHFA, DG and PL.
#'
#' @param chains character vector indicating the FAs to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbTwoChains <- function(chains, lipid) {
  fas <- rep(chains, 2)
  comb <- unique(utils::combn(fas, 2, simplify=F))
  total  <- unique(unlist(lapply(comb, sumChains, n = 2)))
  if (lipid == "PE") {
    fm <- lapply(total, MassPE)
  } else if (lipid == "PG") {
    fm <- lapply(total, MassPG)
  } else if (lipid == "PI") {
    fm <- lapply(total, MassPI)
  } else if (lipid == "PC") {
    fm <- lapply(total, MassPC)
  } else if (lipid == "DG") {
    fm <- lapply(total, MassDG)
  } else if (lipid == "FAHFA") {
    fm <- lapply(total, MassFAHFA)
  } else if (lipid == "PS") {
    fm <- lapply(total, MassPS)
  } else if (lipid == "PA") {
    fm <- lapply(total, MassPA)
  } else if (lipid == "PC-O") {
    fm <- lapply(total, MassPCo)
  } else if (lipid == "PC-P") {
    fm <- lapply(total, MassPCp)
  } else if (lipid == "PE-O") {
    fm <- lapply(total, MassPEo)
  } else if (lipid == "PE-P") {
    fm <- lapply(total, MassPEp)
  } else if (lipid == "MGDG") {
    fm <- lapply(total, MassMGDG)
  } else if (lipid == "DGDG") {
    fm <- lapply(total, MassDGDG)
  } else if (lipid == "BMP") {
    fm <- lapply(total, MassBMP)
  } else if (lipid == "LNAPE") {
    fm <- lapply(total, MassLNAPE)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=total,
                   Mass=as.numeric(as.vector(unlist(lapply(fm, "[[", 2)))),
                   stringsAsFactors = F)
  return(db)
}

# dbPlasmalogens
#' Creation of a database for PCo, PCp, PEo and PEp
#'
#' Creation of a database for PCo, PCp, PEo and PEp
#'
#' @param chains character vector indicating the FAs to be employed for sn1
#' @param chains2 character vector indicating the FAs to be employed for sn2
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbPlasmalogens <- function(chains, chains2, lipid) {
  comb <- vector()
  for (i in chains2) {
    comb <- append(comb, sapply(chains, paste, i, collapse=" "))
  }
  comb <- unique(comb)
  total  <- unique(unlist(lapply(comb, sumChains, n = 2)))
  if (lipid == "PC-O") {
    fm <- lapply(total, MassPCo)
  } else if (lipid == "PC-P") {
    fm <- lapply(total, MassPCp)
  } else if (lipid == "PE-O") {
    fm <- lapply(total, MassPEo)
  } else if (lipid == "PE-P") {
    fm <- lapply(total, MassPEp)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=total,
                   Mass=as.numeric(unlist(lapply(fm, "[[", 2))), stringsAsFactors = F)
  return(db)
}

# dbThreeChains
#' Creation of a database for TG.
#'
#' Creation of a database for TG.
#'
#' @param chains character vector indicating the FAs to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbThreeChains <- function(chains, lipid) {
  fas <- rep(chains, 3)
  comb <- unique(utils::combn(fas, 3, simplify=F))
  total  <- unique(unlist(lapply(comb, sumChains, n = 3)))
  if (lipid == "TG") {
    fm <- lapply(total, MassTG)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=total,
                   Mass=as.numeric(as.vector(unlist(lapply(fm, "[[", 2)))),
                   stringsAsFactors = F)
  return(db)
}

# dbFourChains
#' Creation of a database for C.
#'
#' Creation of a database CL.
#'
#' @param chains character vector indicating the FAs to be employed
#' @param lipid character value indication the class of lipid.
#'
#' @return data frame containing formula, mass and total number of carbons and
#' insaturations
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dbFourChains <- function(chains, lipid) {
  fas <- rep(chains, 4)
  comb <- unique(utils::combn(fas, 4, simplify=F))
  total  <- unique(unlist(lapply(comb, sumChains, n = 4)))
  if (lipid == "CL") {
    fm <- lapply(total, MassCL)
  }
  db <- data.frame(formula=unlist(lapply(fm, "[[", 1)), total=total,
                   Mass=as.numeric(as.vector(unlist(lapply(fm, "[[", 2)))),
                   stringsAsFactors = F)
  return(db)
}


# createLipidDBs
#' Customizable lipid DBs creator
#'
#' It allows to create easy-customizable lipid DBs for annotation with LipidMS
#' package.
#'
#' @param lipid character value indicating the class of lipid. See Details.
#' @param chains character vector indicating the FA chains to be employed
#' @param sph character vector containing the dihydroxy sphingoid bases to be
#' employed if required.
#' @param phyto character vector containing the trihydroxy sphingoid bases to be
#' employed if required.
#' @param alkyl character vector containing the alkyl bases to be employed
#' if required.
#' @param vinyl character vector containing the vinyl ether bases to be employed
#' if required.
#'
#' @return List with the requested dbs (data frames)
#'
#' @details \code{lipidClass} argument needs to be one of the following
#' character values: "Cer_d", "Cer_t", "CerP", "HexCer_d", "HexCer_t",
#' "AHexCer_d", "AHexCer_t", "SHexCer_d", "SHexCer_t", "SM_d", "SM_t", "AcylCer",
#' "Carnitine", "CE", "FA", "HFA", "Sph" (sphingoid bases), "SphP", "MG",
#' "LPA", "LPC", "LPE", "LPG", "LPI", "LPS", "NAGly", "NAGlySer", "NAOrn",
#' "NAE", "FAHFA", "DG", "PC", "PE", "PG", "PI", "PS", "PA", "MGDG",
#' "DGDG", "BMP", "LNAPE", "TG", "CL", "PC-O", "PC-P", "PE-O", "PE-P",
#' "LPA-O", "LPC-O", "LPC-P", "LPE-O", "LPE-P" or "all".
#'
#' @examples
#' \dontrun{
#' fas <- c("8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1",
#' "17:0", "18:0", "18:1", "18:2", "18:3", "18:4", "20:0", "20:1", "20:2",
#' "20:3", "20:4", "20:5", "22:0", "22:1", "22:2", "22:3", "22:4", "22:5",
#' "22:6", "24:0", "24:1", "24:2", "24:3", "24:4","24:5","24:6"," 26:0",
#' "28:0", "30:0", "32:0", "34:0", "36:0")
#'
#' sph <- c("16:0", "16:1", "18:0", "18:1", "18:2", "20:0", "20:1")
#' phyto <- c("18:0")
#' alkyl = c("16:0", "18:0", "18:1")
#' vinyl = c("16:0", "18:0")
#'
#' newdb <- createLipidDBs(lipid = "PC", chains = fas, sph = sph, phyto = phyto)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
createLipidDBs <- function(lipid,
                           chains = c("8:0", "10:0", "12:0", "14:0", "14:1",
                                      "15:0", "16:0", "16:1", "17:0", "18:0",
                                      "18:1", "18:2", "18:3", "18:4", "20:0",
                                      "20:1", "20:2", "20:3", "20:4", "20:5",
                                      "22:0", "22:1", "22:2", "22:3", "22:4",
                                      "22:5", "22:6", "24:0", "24:1", "24:2",
                                      "24:3", "24:4", "24:5", "24:6", "26:0",
                                      "28:0", "30:0", "32:0", "34:0", "36:0"),
                           sph = c("16:0", "16:1", "18:0", "18:1",
                                   "18:2", "20:0", "20:1"),
                           phyto = c("18:0"),
                           alkyl = c("8:0", "10:0", "12:0", "14:0", "14:1",
                                     "15:0", "16:0", "16:1", "17:0", "18:0",
                                     "18:1", "18:2", "18:3", "18:4", "20:0",
                                     "20:1", "20:2", "20:3", "20:4", "20:5",
                                     "22:0", "22:1", "22:2", "22:3", "22:4",
                                     "22:5", "22:6", "24:0", "24:1", "24:2",
                                     "24:3", "24:4", "24:5", "24:6", "26:0",
                                     "28:0", "30:0", "32:0", "34:0", "36:0"),
                           vinyl = c("8:0", "10:0", "12:0", "14:0", "14:1",
                                     "15:0", "16:0", "16:1", "17:0", "18:0",
                                     "18:1", "18:2", "18:3", "18:4", "20:0",
                                     "20:1", "20:2", "20:3", "20:4", "20:5",
                                     "22:0", "22:1", "22:2", "22:3", "22:4",
                                     "22:5", "22:6", "24:0", "24:1", "24:2",
                                     "24:3", "24:4", "24:5", "24:6", "26:0",
                                     "28:0", "30:0", "32:0", "34:0", "36:0")) {

# -------------------------------------------------------------------------

  customizedDataSets <- list()

  if (lipid %in% c("Cer_d", "Cer_t", "CerP_d", "CerP_t", "HexCer_d", "HexCer_t",
                   "AHexCer_d", "AHexCer_t", "SHexCer_d", "SHexCer_t",
                   "SM_d", "SM_t", "AcylCer_d", "AcylCer_t")) {

    if (lipid %in% c("Cer_t", "CerP_t", "HexCer_t", "AHexCer_t", "SHexCer_t",
                     "SM_t", "AcylCer_t")) {
      db <- dbSphingolipids(chains = chains, chains2 = phyto, lipid = lipid)
    } else {
      db <- dbSphingolipids(chains = chains, chains2 = sph, lipid = lipid)
    }

    db <- data.frame(formula = db$formula, total = db$total,
                     Mass = as.numeric(db$Mass),
                     ID = paste(lipid, "(", db$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    customizedDataSets[[lipid]] <- data.frame(formula = db$formula,
                                              total = db$total,
                                              Mass = as.numeric(db$Mass),
                                              ID = db$ID,
                                              stringsAsFactors = FALSE)

  } else if (lipid %in% c("PC-O", "PE-O")) {

    db <- dbPlasmalogens(chains = alkyl, chains2 = chains, lipid = lipid)
    db <- data.frame(formula = db$formula, total = db$total,
                     Mass = as.numeric(db$Mass),
                     ID = paste(lipid, "(", db$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    customizedDataSets[[lipid]] <- data.frame(formula = db$formula,
                                              total = db$total,
                                              Mass = as.numeric(db$Mass),
                                              ID = db$ID,
                                              stringsAsFactors = FALSE)

  } else if (lipid %in% c("PC-P", "PE-P")) {

    db <- dbPlasmalogens(chains = vinyl, chains2 = chains, lipid = lipid)
    db <- data.frame(formula = db$formula, total = db$total,
                     Mass = as.numeric(db$Mass),
                     ID = paste(lipid, "(", db$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    customizedDataSets[[lipid]] <- data.frame(formula = db$formula,
                                              total = db$total,
                                              Mass = as.numeric(db$Mass),
                                              ID = db$ID,
                                              stringsAsFactors = FALSE)

  } else if (lipid %in% c("LPE-O", "LPA-O", "LPC-P", "LPC-O", "LPE-P")) {

    if (lipid %in% c("LPC-O", "LPE-O", "LPA-O")) {
      db <- dbOneChain(chains = alkyl, lipid = lipid)
    } else {
      db <- dbOneChain(chains = vinyl, lipid = lipid)
    }

    db <- data.frame(formula = db$formula, total = db$total,
                     Mass = as.numeric(db$Mass),
                     ID = paste(lipid, "(", db$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    customizedDataSets[[lipid]] <- data.frame(formula = db$formula,
                                              total = db$total,
                                              Mass = as.numeric(db$Mass),
                                              ID = db$ID,
                                              stringsAsFactors = FALSE)

  } else if (lipid %in% c("FA", "HFA", "Carnitine", "LPA", "LPE", "LPG",
                          "LPI", "LPS", "LPC", "MG", "CE",
                          "Sph_d", "Sph_t", "SphP_d", "SphP_t",
                          "NAGly", "NAGlySer", "NAOrn", "NAE")) {

    if (lipid %in% c("Sph_d", "Sph_t", "SphP_d", "SphP_t")) {
      db <- dbOneChain(chains = sph, lipid = lipid)
    } else {
      db <- dbOneChain(chains = chains, lipid = lipid)
    }

    db <- data.frame(formula = db$formula, total = db$total,
                     Mass = as.numeric(db$Mass),
                     ID = paste(lipid, "(", db$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    customizedDataSets[[lipid]] <- data.frame(formula = db$formula,
                                              total = db$total,
                                              Mass = as.numeric(db$Mass),
                                              ID = db$ID,
                                              stringsAsFactors = FALSE)

  } else if (lipid %in% c("DG", "PC", "PE", "PG", "PI", "PS", "FAHFA",
                          "PA", "MGDG", "DGDG", "BMP", "LNAPE")) {

    db <- dbTwoChains(chains = chains, lipid = lipid)
    db <- data.frame(formula = db$formula, total = db$total,
                     Mass = as.numeric(db$Mass),
                     ID = paste(lipid, "(", db$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    customizedDataSets[[lipid]] <- data.frame(formula = db$formula,
                                              total = db$total,
                                              Mass = as.numeric(db$Mass),
                                              ID = db$ID,
                                              stringsAsFactors = FALSE)

  } else if (lipid == "TG") {

    db <- dbThreeChains(chains = chains, lipid = lipid)
    db <- data.frame(formula = db$formula, total = db$total,
                     Mass = as.numeric(db$Mass),
                     ID = paste(lipid, "(", db$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    customizedDataSets[["TG"]] <- data.frame(formula = db$formula,
                                             total = db$total,
                                             Mass = as.numeric(db$Mass),
                                             ID = db$ID,
                                             stringsAsFactors = FALSE)

  } else if (lipid == "CL") {

    db <- dbFourChains(chains = chains, lipid = lipid)
    db <- data.frame(formula = db$formula, total = db$total,
                     Mass = as.numeric(db$Mass),
                     ID = paste(lipid, "(", db$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    customizedDataSets[["CL"]] <- data.frame(formula = db$formula,
                                             total = db$total,
                                             Mass = as.numeric(db$Mass),
                                             ID = db$ID,
                                             stringsAsFactors = FALSE)

  } else if (lipid == "all") {

    cer_d <- dbSphingolipids(chains = chains, chains2 = sph, lipid = "Cer_d")
    cer_d <- data.frame(formula = cer_d$formula, total = cer_d$total,
                        Mass = as.numeric(cer_d$Mass),
                        ID = paste("Cer_d(", cer_d$total, ")", sep = ""),
                        stringsAsFactors = FALSE)

    cer_t <- dbSphingolipids(chains = chains, chains2 = phyto, lipid = "Cer_t")
    cer_t <- data.frame(formula = cer_t$formula, total = cer_t$total,
                        Mass = as.numeric(cer_t$Mass),
                        ID = paste("Cer_t(", cer_t$total, ")", sep = ""),
                        stringsAsFactors = FALSE)

    cerP_d <- dbSphingolipids(chains = chains, chains2 = sph, lipid = "CerP_d")
    cerP_d <- data.frame(formula = cerP_d$formula, total = cerP_d$total,
                       Mass = as.numeric(cerP_d$Mass),
                       ID = paste("CerP_d(", cerP_d$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    cerP_t <- dbSphingolipids(chains = chains, chains2 = phyto, lipid = "CerP_t")
    cerP_t <- data.frame(formula = cerP_t$formula, total = cerP_t$total,
                       Mass = as.numeric(cerP_t$Mass),
                       ID = paste("CerP_t(", cerP_t$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    acylcer_d <- dbSphingolipids(chains = chains, chains2 = sph, lipid = "AcylCer_d")
    acylcer_d <- data.frame(formula = acylcer_d$formula, total = acylcer_d$total,
                          Mass = as.numeric(acylcer_d$Mass),
                          ID = paste("AcylCer_d(", acylcer_d$total, ")", sep = ""),
                          stringsAsFactors = FALSE)

    acylcer_t <- dbSphingolipids(chains = chains, chains2 = phyto, lipid = "AcylCer_t")
    acylcer_t <- data.frame(formula = acylcer_t$formula, total = acylcer_t$total,
                          Mass = as.numeric(acylcer_t$Mass),
                          ID = paste("AcylCer_t(", acylcer_t$total, ")", sep = ""),
                          stringsAsFactors = FALSE)

    hexcer_d <- dbSphingolipids(chains = chains, chains2 = sph, lipid = "HexCer_d")
    hexcer_d <- data.frame(formula = hexcer_d$formula, total = hexcer_d$total,
                           Mass = as.numeric(hexcer_d$Mass),
                           ID = paste("HexCer_d(", hexcer_d$total, ")", sep = ""),
                           stringsAsFactors = FALSE)

    hexcer_t <- dbSphingolipids(chains = chains, chains2 = phyto, lipid = "HexCer_t")
    hexcer_t <- data.frame(formula = hexcer_t$formula, total = hexcer_t$total,
                           Mass = as.numeric(hexcer_t$Mass),
                           ID = paste("HexCer_t(", hexcer_t$total, ")", sep = ""),
                           stringsAsFactors = FALSE)

    ahexcer_d <- dbSphingolipids(chains = chains, chains2 = sph, lipid = "AHexCer_d")
    ahexcer_d <- data.frame(formula = ahexcer_d$formula, total = ahexcer_d$total,
                            Mass = as.numeric(ahexcer_d$Mass),
                            ID = paste("AHexCer_d(", ahexcer_d$total, ")", sep = ""),
                            stringsAsFactors = FALSE)

    ahexcer_t <- dbSphingolipids(chains = chains, chains2 = phyto, lipid = "AHexCer_t")
    ahexcer_t <- data.frame(formula = ahexcer_t$formula, total = ahexcer_t$total,
                            Mass = as.numeric(ahexcer_t$Mass),
                            ID = paste("AHexCer_t(", ahexcer_t$total, ")", sep = ""),
                            stringsAsFactors = FALSE)

    shexcer_d <- dbSphingolipids(chains = chains, chains2 = sph, lipid = "SHexCer_d")
    shexcer_d <- data.frame(formula = shexcer_d$formula, total = shexcer_d$total,
                            Mass = as.numeric(shexcer_d$Mass),
                            ID = paste("SHexCer_d(", shexcer_d$total, ")", sep = ""),
                            stringsAsFactors = FALSE)

    shexcer_t <- dbSphingolipids(chains = chains, chains2 = phyto, lipid = "SHexCer_t")
    shexcer_t <- data.frame(formula = shexcer_t$formula, total = shexcer_t$total,
                            Mass = as.numeric(shexcer_t$Mass),
                            ID = paste("SHexCer_t(", shexcer_t$total, ")", sep = ""),
                            stringsAsFactors = FALSE)

    sm_d <- dbSphingolipids(chains = chains, chains2 = sph, lipid = "SM_d")
    sm_d <- data.frame(formula = sm_d$formula, total = sm_d$total,
                     Mass = as.numeric(sm_d$Mass),
                     ID = paste("SM_d(", sm_d$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    sm_t <- dbSphingolipids(chains = chains, chains2 = phyto, lipid = "SM_t")
    sm_t <- data.frame(formula = sm_t$formula, total = sm_t$total,
                       Mass = as.numeric(sm_t$Mass),
                       ID = paste("SM_t(", sm_t$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    fa <- dbOneChain(chains = chains, lipid = "FA")
    fa <- data.frame(formula = fa$formula, total = fa$total,
                     Mass = as.numeric(fa$Mass),
                     ID = paste("FA(", fa$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    hfa <- dbOneChain(chains = chains, lipid = "HFA")
    hfa <- data.frame(formula = hfa$formula, total = hfa$total,
                      Mass = as.numeric(hfa$Mass),
                      ID = paste("HFA(", hfa$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    carnitine <- dbOneChain(chains = chains, lipid = "Carnitine")
    carnitine <- data.frame(formula = carnitine$formula, total = carnitine$total,
                            Mass = as.numeric(carnitine$Mass),
                            ID = paste("Carnitine(", carnitine$total, ")", sep = ""),
                            stringsAsFactors = FALSE)

    CE <- dbOneChain(chains = chains, lipid = "CE")
    CE <- data.frame(formula = CE$formula, total = CE$total,
                     Mass = as.numeric(CE$Mass),
                     ID = paste("CE(", CE$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    mg <- dbOneChain(chains = chains, lipid = "MG")
    mg <- data.frame(formula = mg$formula, total = mg$total,
                     Mass = as.numeric(mg$Mass),
                     ID = paste("MG(", mg$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    sphdb_d <- dbOneChain(chains = sph, lipid = "Sph_d")
    sphdb_d <- data.frame(formula = sphdb_d$formula, total = sphdb_d$total,
                        Mass = as.numeric(sphdb_d$Mass),
                        ID = paste("Sph_d(", sphdb_d$total, ")", sep = ""),
                        stringsAsFactors = FALSE)

    sphdb_t <- dbOneChain(chains = sph, lipid = "Sph_t")
    sphdb_t <- data.frame(formula = sphdb_t$formula, total = sphdb_t$total,
                        Mass = as.numeric(sphdb_t$Mass),
                        ID = paste("Sph_t(", sphdb_t$total, ")", sep = ""),
                        stringsAsFactors = FALSE)

    nlsph <- sphdb_d[, 1:3]
    nlsph$Mass <- nlsph$Mass - 44.05

    sphP_d <- dbOneChain(chains = sph, lipid = "SphP_d")
    sphP_d <- data.frame(formula = sphP_d$formula, total = sphP_d$total,
                       Mass = as.numeric(sphP_d$Mass),
                       ID = paste("SphP_d(", sphP_d$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    sphP_t <- dbOneChain(chains = sph, lipid = "SphP_t")
    sphP_t <- data.frame(formula = sphP_t$formula, total = sphP_t$total,
                       Mass = as.numeric(sphP_t$Mass),
                       ID = paste("SphP_t(", sphP_t$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    lpa <- dbOneChain(chains = chains, lipid = "LPA")
    lpa <- data.frame(formula = lpa$formula, total = lpa$total,
                      Mass = as.numeric(lpa$Mass),
                      ID = paste("LPA(", lpa$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lpe <- dbOneChain(chains = chains, lipid = "LPE")
    lpe <- data.frame(formula = lpe$formula, total = lpe$total,
                      Mass = as.numeric(lpe$Mass),
                      ID = paste("LPE(", lpe$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lpg <- dbOneChain(chains = chains, lipid = "LPG")
    lpg <- data.frame(formula = lpg$formula, total = lpg$total,
                      Mass = as.numeric(lpg$Mass),
                      ID = paste("LPG(", lpg$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lpi <- dbOneChain(chains = chains, lipid = "LPI")
    lpi <- data.frame(formula = lpi$formula, total = lpi$total,
                      Mass = as.numeric(lpi$Mass),
                      ID = paste("LPI(", lpi$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lps <- dbOneChain(chains = chains, lipid = "LPS")
    lps <- data.frame(formula = lps$formula, total = lps$total,
                      Mass = as.numeric(lps$Mass),
                      ID = paste("LPS(", lps$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lpc <- dbOneChain(chains = chains, lipid = "LPC")
    lpc <- data.frame(formula = lpc$formula, total = lpc$total,
                      Mass = as.numeric(lpc$Mass),
                      ID = paste("LPC(", lpc$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    nagly <- dbOneChain(chains = chains, lipid = "NAGly")
    nagly <- data.frame(formula = nagly$formula, total = nagly$total,
                        Mass = as.numeric(nagly$Mass),
                        ID = paste("NAGly(", nagly$total, ")", sep = ""),
                        stringsAsFactors = FALSE)

    naglyser <- dbOneChain(chains = chains, lipid = "NAGlySer")
    naglyser <- data.frame(formula = naglyser$formula, total = naglyser$total,
                           Mass = as.numeric(naglyser$Mass),
                           ID = paste("NAGlySer(", naglyser$total, ")", sep = ""),
                           stringsAsFactors = FALSE)

    naorn <- dbOneChain(chains = chains, lipid = "NAOrn")
    naorn <- data.frame(formula = naorn$formula, total = naorn$total,
                        Mass = as.numeric(naorn$Mass),
                        ID = paste("NAOrn(", naorn$total, ")", sep = ""),
                        stringsAsFactors = FALSE)

    nae <- dbOneChain(chains = chains, lipid = "NAE")
    nae <- data.frame(formula = nae$formula, total = nae$total,
                      Mass = as.numeric(nae$Mass),
                      ID = paste("NAE(", nae$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    pc <- dbTwoChains(chains = chains, lipid = "PC")
    pc <- data.frame(formula = pc$formula, total = pc$total,
                     Mass = as.numeric(pc$Mass),
                     ID = paste("PC(", pc$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    pe <- dbTwoChains(chains = chains, lipid = "PE")
    pe <- data.frame(formula = pe$formula, total = pe$total,
                     Mass = as.numeric(pe$Mass),
                     ID = paste("PE(", pe$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    pg <- dbTwoChains(chains = chains, lipid = "PG")
    pg <- data.frame(formula = pg$formula, total = pg$total,
                     Mass = as.numeric(pg$Mass),
                     ID = paste("PG(", pg$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    pi <- dbTwoChains(chains = chains, lipid = "PI")
    pi <- data.frame(formula = pi$formula, total = pi$total,
                     Mass = as.numeric(pi$Mass),
                     ID = paste("PI(", pi$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    ps <- dbTwoChains(chains = chains, lipid = "PS")
    ps <- data.frame(formula = ps$formula, total = ps$total,
                     Mass = as.numeric(ps$Mass),
                     ID = paste("PS(", ps$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    pa <- dbTwoChains(chains = chains, lipid = "PA")
    pa <- data.frame(formula = pa$formula, total = pa$total,
                     Mass = as.numeric(pa$Mass),
                     ID = paste("PA(", pa$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    fahfa <- dbTwoChains(chains = chains, lipid = "FAHFA")
    fahfa <- data.frame(formula = fahfa$formula, total = fahfa$total,
                        Mass = as.numeric(fahfa$Mass),
                        ID = paste("FAHFA(", fahfa$total, ")", sep = ""),
                        stringsAsFactors = FALSE)

    dg <- dbTwoChains(chains = chains, lipid = "DG")
    dg <- data.frame(formula = dg$formula, total = dg$total,
                     Mass = as.numeric(dg$Mass),
                     ID = paste("DG(", dg$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    mgdg <- dbTwoChains(chains = chains, lipid = "MGDG")
    mgdg <- data.frame(formula = mgdg$formula, total = mgdg$total,
                       Mass = as.numeric(mgdg$Mass),
                       ID = paste("MGDG(", mgdg$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    dgdg <- dbTwoChains(chains = chains, lipid = "DGDG")
    dgdg <- data.frame(formula = dgdg$formula, total = dgdg$total,
                       Mass = as.numeric(dgdg$Mass),
                       ID = paste("DGDG(", dgdg$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    bmp <- dbTwoChains(chains = chains, lipid = "BMP")
    bmp <- data.frame(formula = bmp$formula, total = bmp$total,
                      Mass = as.numeric(bmp$Mass),
                      ID = paste("BMP(", bmp$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lnape <- dbTwoChains(chains = chains, lipid = "LNAPE")
    lnape <- data.frame(formula = lnape$formula, total = lnape$total,
                        Mass = as.numeric(lnape$Mass),
                        ID = paste("LNAPE(", lnape$total, ")", sep = ""),
                        stringsAsFactors = FALSE)

    tg <- dbThreeChains(chains = chains, lipid = "TG")
    tg <- data.frame(formula = tg$formula, total = tg$total,
                     Mass = as.numeric(tg$Mass),
                     ID = paste("TG(", tg$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    cl <- dbFourChains(chains = chains, lipid = "CL")
    cl <- data.frame(formula = cl$formula, total = cl$total,
                     Mass = as.numeric(cl$Mass),
                     ID = paste("CL(", cl$total, ")", sep = ""),
                     stringsAsFactors = FALSE)

    peo <- dbPlasmalogens(chains = alkyl, chains2 = chains, lipid = "PE-O")
    peo <- data.frame(formula = peo$formula, total = peo$total,
                      Mass = as.numeric(peo$Mass),
                      ID = paste("PE-O(", peo$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lpeo <- dbOneChain(chains = alkyl, lipid = "LPE-O")
    lpeo <- data.frame(formula = lpeo$formula, total = lpeo$total,
                       Mass = as.numeric(lpeo$Mass),
                       ID = paste("LPE-O(", lpeo$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    lpao <- dbOneChain(chains = alkyl, lipid = "LPA-O")
    lpao <- data.frame(formula = lpao$formula, total = lpao$total,
                       Mass = as.numeric(lpao$Mass),
                       ID = paste("LPA-O(", lpao$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    pco <- dbPlasmalogens(chains = alkyl, chains2 = chains, lipid = "PC-O")
    pco <- data.frame(formula = pco$formula, total = pco$total,
                      Mass = as.numeric(pco$Mass),
                      ID = paste("PC-O(", pco$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lpco <- dbOneChain(chains = alkyl, lipid = "LPC-O")
    lpco <- data.frame(formula = lpco$formula, total = lpco$total,
                       Mass = as.numeric(lpco$Mass),
                       ID = paste("LPC-O(", lpco$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    pep <- dbPlasmalogens(chains = vinyl, chains2 = chains, lipid = "PE-P")
    pep <- data.frame(formula = pep$formula, total = pep$total,
                      Mass = as.numeric(pep$Mass),
                      ID = paste("PE-P(", pep$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lpep <- dbOneChain(chains = vinyl, lipid = "LPE-P")
    lpep <- data.frame(formula = lpep$formula, total = lpep$total,
                       Mass = as.numeric(lpep$Mass),
                       ID = paste("LPE-P(", lpep$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    pcp <- dbPlasmalogens(chains = vinyl, chains2 = chains, lipid = "PC-P")
    pcp <- data.frame(formula = pcp$formula, total = pcp$total,
                      Mass = as.numeric(pcp$Mass),
                      ID = paste("PC-P(", pcp$total, ")", sep = ""),
                      stringsAsFactors = FALSE)

    lpcp <- dbOneChain(chains = vinyl, lipid = "LPC-P")
    lpcp <- data.frame(formula = lpcp$formula, total = lpcp$total,
                       Mass = as.numeric(lpcp$Mass),
                       ID = paste("LPC-P(", lpcp$total, ")", sep = ""),
                       stringsAsFactors = FALSE)

    customizedDataSets[["Cer_d"]] <- cer_d
    customizedDataSets[["Cer_t"]] <- cer_t
    customizedDataSets[["CerP_d"]] <- cerP_d
    customizedDataSets[["CerP_t"]] <- cerP_t
    customizedDataSets[["AcylCer_d"]] <- acylcer_d
    customizedDataSets[["AcylCer_t"]] <- acylcer_t
    customizedDataSets[["HexCer_d"]] <- hexcer_d
    customizedDataSets[["HexCer_t"]] <- hexcer_t
    customizedDataSets[["AHexCer_d"]] <- ahexcer_d
    customizedDataSets[["AHexCer_t"]] <- ahexcer_t
    customizedDataSets[["SHexCer_d"]] <- shexcer_d
    customizedDataSets[["SHexCer_t"]] <- shexcer_t
    customizedDataSets[["SM_d"]] <- sm_d
    customizedDataSets[["SM_t"]] <- sm_t
    customizedDataSets[["FA"]] <- fa
    customizedDataSets[["HFA"]] <- hfa
    customizedDataSets[["Carnitine"]] <- carnitine
    customizedDataSets[["LPA"]] <- lpa
    customizedDataSets[["LPE"]] <- lpe
    customizedDataSets[["LPG"]] <- lpg
    customizedDataSets[["LPI"]] <- lpi
    customizedDataSets[["LPS"]] <- lps
    customizedDataSets[["LPC"]] <- lpc
    customizedDataSets[["MG"]] <- mg
    customizedDataSets[["CE"]] <- CE
    customizedDataSets[["Sph_d"]] <- sphdb_d
    customizedDataSets[["Sph_t"]] <- sphdb_t
    customizedDataSets[["SphP_d"]] <- sphP_d
    customizedDataSets[["SphP_t"]] <- sphP_t
    customizedDataSets[["NAGly"]] <- nagly
    customizedDataSets[["NAGlySer"]] <- naglyser
    customizedDataSets[["NAOrn"]] <- naorn
    customizedDataSets[["NAE"]] <- nae
    customizedDataSets[["FAHFA"]] <- fahfa
    customizedDataSets[["PE"]] <- pe
    customizedDataSets[["PG"]] <- pg
    customizedDataSets[["PI"]] <- pi
    customizedDataSets[["PS"]] <- ps
    customizedDataSets[["PA"]] <- pa
    customizedDataSets[["PC"]] <- pc
    customizedDataSets[["DG"]] <- dg
    customizedDataSets[["MGDG"]] <- mgdg
    customizedDataSets[["DGDG"]] <- dgdg
    customizedDataSets[["BMP"]] <- bmp
    customizedDataSets[["LNAPE"]] <- lnape
    customizedDataSets[["TG"]] <- tg
    customizedDataSets[["CL"]] <- cl
    customizedDataSets[["PE-O"]] <- peo
    customizedDataSets[["LPE-O"]] <- lpeo
    customizedDataSets[["LPA-O"]] <- lpao
    customizedDataSets[["PC-O"]] <- pco
    customizedDataSets[["LPC-O"]] <- lpco
    customizedDataSets[["PE-P"]] <- pep
    customizedDataSets[["LPE-P"]] <- lpep
    customizedDataSets[["PC-P"]] <- pcp
    customizedDataSets[["LPC-P"]] <- lpcp
  }

  return(customizedDataSets)
}

# assignDBs
#' Load LipidRT default data bases
#'
#' Load all LipidRT default data bases required to run identification functions.
#'
#' @return List of data frames.
#'
#' @examples
#' \dontrun{
#' dbs <- assignDBs()
#' }
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
assignDBs <- function() {
  dbs <- list()

  dbs[["Cer_d"]] <- LipidRT::cerddb
  dbs[["Cer_t"]] <- LipidRT::certdb
  dbs[["CerP_d"]] <- LipidRT::cerPddb
  dbs[["CerP_t"]] <- LipidRT::cerPtdb
  dbs[["AcylCer_d"]] <- LipidRT::acylcerddb
  dbs[["AcylCer_t"]] <- LipidRT::acylcertdb
  dbs[["HexCer_d"]] <- LipidRT::hexcerddb
  dbs[["HexCer_t"]] <- LipidRT::hexcertdb
  dbs[["AHexCer_d"]] <- LipidRT::ahexcerddb
  dbs[["AHexCer_t"]] <- LipidRT::ahexcertdb
  dbs[["SHexCer_d"]] <- LipidRT::shexcerddb
  dbs[["SHexCer_t"]] <- LipidRT::shexcertdb
  dbs[["SM_d"]] <- LipidRT::smddb
  dbs[["SM_t"]] <- LipidRT::smtdb

  dbs[["FA"]] <- LipidRT::fadb
  dbs[["HFA"]] <- LipidRT::hfadb
  dbs[["Carnitine"]] <- LipidRT::carnitinesdb
  dbs[["CE"]] <- LipidRT::CEdb
  dbs[["Sph_d"]] <- LipidRT::sphddb
  dbs[["Sph_t"]] <- LipidRT::sphtdb
  dbs[["SphP_d"]] <- LipidRT::sphPddb
  dbs[["SphP_t"]] <- LipidRT::sphPtdb

  dbs[["LPA"]] <- LipidRT::lysopadb
  dbs[["LPA-O"]] <- LipidRT::lysopaodb
  dbs[["LPE"]] <- LipidRT::lysopedb
  dbs[["LPE-O"]] <- LipidRT::lysopeodb
  dbs[["LPE-P"]] <- LipidRT::lysopepdb
  dbs[["LPG"]] <- LipidRT::lysopgdb
  dbs[["LPI"]] <- LipidRT::lysopidb
  dbs[["LPS"]] <- LipidRT::lysopsdb
  dbs[["LPC"]] <- LipidRT::lysopcdb
  dbs[["LPC-O"]] <- LipidRT::lysopcodb
  dbs[["LPC-P"]] <- LipidRT::lysopcpdb

  dbs[["MG"]] <- LipidRT::mgdb
  dbs[["DG"]] <- LipidRT::dgdb
  dbs[["TG"]] <- LipidRT::tgdb
  dbs[["CL"]] <- LipidRT::cldb

  dbs[["FAHFA"]] <- LipidRT::fahfadb
  dbs[["NAGly"]] <- LipidRT::naglydb
  dbs[["NAGlySer"]] <- LipidRT::naglyserdb
  dbs[["NAOrn"]] <- LipidRT::naorndb
  dbs[["NAE"]] <- LipidRT::naedb

  dbs[["PC"]] <- LipidRT::pcdb
  dbs[["PC-O"]] <- LipidRT::pcodb
  dbs[["PC-P"]] <- LipidRT::pcpdb
  dbs[["PE"]] <- LipidRT::pedb
  dbs[["PE-O"]] <- LipidRT::peodb
  dbs[["PE-P"]] <- LipidRT::pepdb
  dbs[["PG"]] <- LipidRT::pgdb
  dbs[["PI"]] <- LipidRT::pidb
  dbs[["PS"]] <- LipidRT::psdb
  dbs[["PA"]] <- LipidRT::padb
  dbs[["MGDG"]] <- LipidRT::mgdgdb
  dbs[["DGDG"]] <- LipidRT::dgdgdb
  dbs[["BMP"]] <- LipidRT::bmpdb
  dbs[["LNAPE"]] <- LipidRT::lnapedb

  dbs[["adductsTable"]] <- LipidRT::adductsTable

  return(dbs)
}

# defaultfiledbs
#' Create a default editable file for modifying lipid classes, adducts, and chains
#'
#' This function generates a template data frame (\code{filedbs}) containing
#' default lipid classes, associated adducts, and fatty acyl chains. The
#' resulting file can be used to customize lipid annotations in downstream
#' RT-based analyses.
#'
#' @param sep Character. Column delimiter for the output file (default: ",").
#' @param path Character. Optional path to save the resulting file. If
#' \code{NULL}, the data frame is returned without saving.
#'
#' @return \code{filedbs} data frame
#'
#' @examples
#' \dontrun{
#' filedbs <- defaultfiledbs(sep = ",")
#' }
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
defaultfiledbs <- function(sep = ",", path = NULL) {

  # default chains
  fas = c("8:0", "10:0", "12:0", "14:0", "14:1",
             "15:0", "16:0", "16:1", "17:0", "18:0",
             "18:1", "18:2", "18:3", "18:4", "20:0",
             "20:1", "20:2", "20:3", "20:4", "20:5",
             "22:0", "22:1", "22:2", "22:3", "22:4",
             "22:5", "22:6", "24:0", "24:1", "24:2",
             "24:3", "24:4", "24:5", "24:6", "26:0",
             "28:0", "30:0", "32:0", "34:0", "36:0")
  sph = c("16:0", "16:1", "18:0", "18:1",
          "18:2", "20:0", "20:1")
  phyto = c("18:0")
  alkyl = fas
  vinyl = fas

  lipidsubclasses <- list(
    "Cer_d" = list(dbs = "cerddb",
                   adducts_pos = "M+H-H2O; M+Na; M+H",
                   adducts_neg = "M-H; M+CH3COO; M+CH3COO-CH3",
                   fas = fas,
                   sph = sph,
                   alkyl = NA,
                   vinyl = NA),

    "Cer_t" = list(dbs = "certdb",
                   adducts_pos = "M+H-H2O; M+Na; M+H",
                   adducts_neg = "M-H; M+CH3COO; M+CH3COO-CH3",
                   fas = fas,
                   sph = phyto,
                   alkyl = NA,
                   vinyl = NA),

    "CerP_d" = list(dbs = "cerPddb",
                    adducts_pos = "M+H",
                    adducts_neg = "M-H",
                    fas = fas,
                    sph = sph,
                    alkyl = NA,
                    vinyl = NA),

    "CerP_t" = list(dbs = "cerPtdb",
                    adducts_pos = "M+H",
                    adducts_neg = "M-H",
                    fas = fas,
                    sph = phyto,
                    alkyl = NA,
                    vinyl = NA),

    "AcylCer_d" = list(dbs = "acylcerddb",
                       adducts_pos = "M+H; M+H-H2O; M+Na",
                       adducts_neg = "M-H; M+CH3COO; M+CH3COO-CH3",
                       fas = fas,
                       sph = sph,
                       alkyl = NA,
                       vinyl = NA),

    "AcylCer_t" = list(dbs = "acylcertdb",
                       adducts_pos = "M+H; M+H-H2O; M+Na",
                       adducts_neg = "M-H; M+CH3COO; M+CH3COO-CH3",
                       fas = fas,
                       sph = phyto,
                       alkyl = NA,
                       vinyl = NA),

    "HexCer_d" = list(dbs = "hexcerddb",
                      adducts_pos = "M+H; M+Na",
                      adducts_neg = "M-H; M+CH3COO; M+CH3COO-CH3",
                      fas = fas,
                      sph = sph,
                      alkyl = NA,
                      vinyl = NA),

    "HexCer_t" = list(dbs = "hexcertdb",
                      adducts_pos = "M+H; M+Na",
                      adducts_neg = "M-H; M+CH3COO; M+CH3COO-CH3",
                      fas = fas,
                      sph = phyto,
                      alkyl = NA,
                      vinyl = NA),

    "AHexCer_d" = list(dbs = "ahexcerddb",
                       adducts_pos = "M+H; M+Na",
                       adducts_neg = "M-H; M+CH3COO; M+CH3COO-CH3",
                       fas = fas,
                       sph = sph,
                       alkyl = NA,
                       vinyl = NA),

    "AHexCer_t" = list(dbs = "ahexcertdb",
                       adducts_pos = "M+H; M+Na",
                       adducts_neg = "M-H; M+CH3COO; M+CH3COO-CH3",
                       fas = fas,
                       sph = phyto,
                       alkyl = NA,
                       vinyl = NA),

    "SHexCer_d" = list(dbs = "shexcerddb",
                       adducts_pos = "M+H; M+Na",
                       adducts_neg = "M-H",
                       fas = fas,
                       sph = sph,
                       alkyl = NA,
                       vinyl = NA),

    "SHexCer_t" = list(dbs = "shexcertdb",
                       adducts_pos = "M+H; M+Na",
                       adducts_neg = "M-H",
                       fas = fas,
                       sph = phyto,
                       alkyl = NA,
                       vinyl = NA),

    "SM_d" = list(dbs = "smddb",
                adducts_pos = "M+H; M+Na",
                adducts_neg = "M+CH3COO; M-CH3; M+CH3COO-CH3",
                fas = fas,
                sph = sph,
                alkyl = NA,
                vinyl = NA),

    "SM_t" = list(dbs = "smtdb",
                adducts_pos = "M+H; M+Na",
                adducts_neg = "M+CH3COO; M-CH3; M+CH3COO-CH3",
                fas = fas,
                sph = sph,
                alkyl = NA,
                vinyl = NA),

    "FA" = list(dbs = "fadb",
                adducts_pos = NA,
                adducts_neg = "M-H; 2M-H",
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "HFA" = list(dbs = "hfadb",
                 adducts_pos = NA,
                 adducts_neg = "M-H",
                 fas = fas,
                 sph = NA,
                 alkyl = NA,
                 vinyl = NA),

    "Carnitine" = list(dbs = "carnitinesdb",
                       adducts_pos = "M+H; M+Na",
                       adducts_neg = NA,
                       fas = fas,
                       sph = NA,
                       alkyl = NA,
                       vinyl = NA),

    "LPA" = list(dbs = "lysopadb",
                 adducts_pos = NA,
                 adducts_neg = "M-H",
                 fas = fas,
                 sph = NA,
                 alkyl = NA,
                 vinyl = NA),

    "LPE" = list(dbs = "lysopedb",
                 adducts_pos = "M+H; M+Na",
                 adducts_neg = "M-H",
                 fas = fas,
                 sph = NA,
                 alkyl = NA,
                 vinyl = NA),

    "LPG" = list(dbs = "lysopgdb",
                 adducts_pos = NA,
                 adducts_neg = "M-H",
                 fas = fas,
                 sph = NA,
                 alkyl = NA,
                 vinyl = NA),

    "LPI" = list(dbs = "lysopidb",
                 adducts_pos = NA,
                 adducts_neg = "M-H",
                 fas = fas,
                 sph = NA,
                 alkyl = NA,
                 vinyl = NA),

    "LPS" = list(dbs = "lysopsdb",
                 adducts_pos = NA,
                 adducts_neg = "M-H; M+Na-2H",
                 fas = fas,
                 sph = NA,
                 alkyl = NA,
                 vinyl = NA),

    "LPC" = list(dbs = "lysopcdb",
                 adducts_pos = "M+H; M+Na",
                 adducts_neg = "M+CH3COO; M-CH3; M+CH3COO-CH3",
                 fas = fas,
                 sph = NA,
                 alkyl = NA,
                 vinyl = NA),

    "MG" = list(dbs = "mgdb",
                adducts_pos = "M+H-H2O; M+NH4; M+Na",
                adducts_neg = NA,
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "CE" = list(dbs = "CEdb",
                adducts_pos = "2M+NH4; 2M+Na; M+NH4; M+Na",
                adducts_neg = NA,
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "Sph_d" = list(dbs = "sphddb",
                 adducts_pos = "M+H",
                 adducts_neg = "M-H",
                 fas = NA,
                 sph = sph,
                 alkyl = NA,
                 vinyl = NA),

    "Sph_t" = list(dbs = "sphtpb",
                   adducts_pos = "M+H",
                   adducts_neg = "M-H",
                   fas = NA,
                   sph = sph,
                   alkyl = NA,
                   vinyl = NA),

    "SphP_d" = list(dbs = "sphPddb",
                    adducts_pos = "M+H",
                    adducts_neg = "M-H",
                    fas = NA,
                    sph = sph,
                    alkyl = NA,
                    vinyl = NA),

    "SphP_t" = list(dbs = "sphPtdb",
                    adducts_pos = "M+H",
                    adducts_neg = "M-H",
                    fas = NA,
                    sph = sph,
                    alkyl = NA,
                    vinyl = NA),

    "NAGly" = list(dbs = "naglydb",
                   adducts_pos = "M+H; M+Na",
                   adducts_neg = "M-H",
                   fas = fas,
                   sph = NA,
                   alkyl = NA,
                   vinyl = NA),

    "NAGlySer" = list(dbs = "naglyserdb",
                      adducts_pos = "M+H; M+Na",
                      adducts_neg = "M-H",
                      fas = fas,
                      sph = NA,
                      alkyl = NA,
                      vinyl = NA),

    "NAOrn" = list(dbs = "naorndb",
                   adducts_pos = "M+H; M+Na",
                   adducts_neg = "M-H",
                   fas = fas,
                   sph = NA,
                   alkyl = NA,
                   vinyl = NA),

    "NAE" = list(dbs = "naedb",
                 adducts_pos = "M+H; M+Na",
                 adducts_neg = "M-H",
                 fas = fas,
                 sph = NA,
                 alkyl = NA,
                 vinyl = NA),

    "DG" = list(dbs = "dgdb",
                adducts_pos = "M+H-H2O; M+NH4; M+Na",
                adducts_neg = NA,
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "PC" = list(dbs = "pcdb",
                adducts_pos = "M+H; M+Na",
                adducts_neg = "M+CH3COO; M-CH3; M+CH3COO-CH3",
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "PE" = list(dbs = "pedb",
                adducts_pos = "M+H; M+Na",
                adducts_neg = "M-H",
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "PG" = list(dbs = "pgdb",
                adducts_pos = "M+H; M+NH4; M+Na",
                adducts_neg = "M-H",
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "PS" = list(dbs = "psdb",
                adducts_pos = NA,
                adducts_neg = "M-H; M+Na-2H",
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "PI" = list(dbs = "pidb",
                adducts_pos = "M+H; M+NH4; M+Na",
                adducts_neg = "M-H",
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "PA" = list(dbs = "padb",
                adducts_pos = NA,
                adducts_neg = "M-H",
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),


    "FAHFA" = list(dbs = "fahfadb",
                   adducts_pos = NA,
                   adducts_neg = "M-H",
                   fas = fas,
                   sph = NA,
                   alkyl = NA,
                   vinyl = NA),

    "MGDG" = list(dbs = "mgdgdb",
                  adducts_pos = "M+NH4; M+Na",
                  adducts_neg = "M-H",
                  fas = fas,
                  sph = NA,
                  alkyl = NA,
                  vinyl = NA),

    "DGDG" = list(dbs = "dgdgdb",
                  adducts_pos = "M+NH4; M+Na",
                  adducts_neg = "M-H",
                  fas = fas,
                  sph = NA,
                  alkyl = NA,
                  vinyl = NA),

    "BMP" = list(dbs = "bmpdb",
                 adducts_pos = "M+NH4; M+Na",
                 adducts_neg = "M-H",
                 fas = fas,
                 sph = NA,
                 alkyl = NA,
                 vinyl = NA),

    "LNAPE" = list(dbs = "lnapedb",
                   adducts_pos = "M+H; M+Na",
                   adducts_neg = "M-H",
                   fas = fas,
                   sph = NA,
                   alkyl = NA,
                   vinyl = NA),

    "TG" = list(dbs = "tgdb",
                adducts_pos = "M+NH4; M+Na",
                adducts_neg = NA,
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "CL" = list(dbs = "cldb",
                adducts_pos = NA,
                adducts_neg = "M-H; M+Na-2H",
                fas = fas,
                sph = NA,
                alkyl = NA,
                vinyl = NA),

    "PE-O" = list(dbs = "peodb",
                  adducts_pos = "M+H; M+Na",
                  adducts_neg = "M-H; M+NaCH3COO",
                  fas = fas,
                  sph = NA,
                  alkyl = alkyl,
                  vinyl = NA),

    "PC-O" = list(dbs = "pcodb",
                  adducts_pos = "M+H; M+Na",
                  adducts_neg = "M+CH3COO; M-CH3; M+CH3COO-CH3",
                  fas = fas,
                  sph = NA,
                  alkyl = alkyl,
                  vinyl = NA),

    "LPE-O" = list(dbs = "lysopeodb",
                   adducts_pos = "M+H; M+Na",
                   adducts_neg = "M-H; M+NaCH3COO",
                   fas = NA,
                   sph = NA,
                   alkyl = alkyl,
                   vinyl = NA),

    "LPC-O" = list(dbs = "lysopcodb",
                   adducts_pos = "M+H; M+Na",
                   adducts_neg = "M+CH3COO; M-CH3; M+CH3COO-CH3",
                   fas = NA,
                   sph = NA,
                   alkyl = alkyl,
                   vinyl = NA),

    "LPA-O" = list(dbs = "lysopaodb",
                   adducts_pos = NA,
                   adducts_neg = "M-H",
                   fas = NA,
                   sph = NA,
                   alkyl = alkyl,
                   vinyl = NA),

    "PE-P" = list(dbs = "pepdb",
                  adducts_pos = "M+H; M+Na",
                  adducts_neg = "M-H; M+NaCH3COO",
                  fas = fas,
                  sph = NA,
                  alkyl = NA,
                  vinyl = vinyl),

    "PC-P" = list(dbs = "pcpdb",
                  adducts_pos = "M+H; M+Na",
                  adducts_neg = "M+CH3COO; M-CH3; M+CH3COO-CH3",
                  fas = fas,
                  sph = NA,
                  alkyl = NA,
                  vinyl = vinyl),

    "LPE-P" = list(dbs = "lysopepdb",
                   adducts_pos = "M+H; M+Na",
                   adducts_neg = "M-H; M+NaCH3COO",
                   fas = NA,
                   sph = NA,
                   alkyl = NA,
                   vinyl = vinyl),

    "LPC-P" = list(dbs = "lysopcpdb",
                   adducts_pos = "M+H; M+Na",
                   adducts_neg = "M+CH3COO; M-CH3; M+CH3COO-CH3",
                   fas = NA,
                   sph = NA,
                   alkyl = NA,
                   vinyl = vinyl)
  )

  df <- do.call(rbind, lapply(names(lipidsubclasses), function(subclass_name) {

    class_data <- lipidsubclasses[[subclass_name]]

    data.frame(
      classes = subclass_name,
      adducts_pos = safe_paste(class_data$adducts_pos),
      adducts_neg = safe_paste(class_data$adducts_neg),
      fas = safe_paste(class_data$fas),
      sph = safe_paste(class_data$sph),
      alkyl = safe_paste(class_data$alkyl),
      vinyl = safe_paste(class_data$vinyl),
      stringsAsFactors = FALSE
    )
  }))

  df[df == "NA" | is.na(df)] <- ""
  df$adducts_pos <- gsub("; ", "\\|", df$adducts_pos)
  df$adducts_neg <- gsub("; ", "\\|", df$adducts_neg)
  df <- df[df$adducts_pos != "" | df$adducts_neg != "", ]

  df$fas <- gsub("; ", "\\|", df$fas)
  df$sph <- gsub("; ", "\\|", df$sph)
  df$alkyl <- gsub("; ", "\\|", df$alkyl)
  df$vinyl <- gsub("; ", "\\|", df$vinyl)


  if (!is.null(path)) {
    if (is.character(path)) {
      if (grepl("\\.csv$", path)) {
        write.csv(df, file = path, row.names = FALSE)
      } else {
        if (!grepl("/$", path)) {
          path <- paste0(path, "/")
        }
        write.table(df, file = paste0(path, "filedbs.csv"), sep = sep,
                    row.names = FALSE, col.names = TRUE, quote = TRUE)
      }
    } else {
      stop("The 'path' argument must be a character string.")
    }
  }

  return(df)
}

safe_paste <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  paste(x, collapse = "|")
}
