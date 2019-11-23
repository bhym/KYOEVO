#DUMMY LIST OF INSTRUCTIONS
DummyList <- list(
DummyInstruction0 = function(x) {return(x)},
DummyInstruction1 = function(x) {return(paste(x, "b", sep=""))},
DummyInstruction2 = function(x) {return(paste(x, "a", sep=""))},
DummyInstruction3 = function(x) {return(paste(x, "o", sep=""))}
)

#THE LIST CONTAINING WRAPPERS TO THE FUNDAMENTAL INSTRUCTIONS
InstructionsList <- lapply(DummyList, function(x) {
                                        Instruction_Class$new(instruction=x)
                                      }
                          )
#names <- paste ("Fun",1:3, sep="")
#invisible(mapply(function(x,y) { x$name <- y}, InstructionsList,names))
names(InstructionsList) = paste("Instruction", seq_along(DummyList), sep="")

#DEFINITION OF THE ALTERNATIVES FOR EACH INSTRUCTION
InstructionsList[[1]]$set_alt(InstructionsList[-1])
InstructionsList[[2]]$set_alt(InstructionsList[-2])
InstructionsList[[3]]$set_alt(InstructionsList[-3])
InstructionsList[[4]]$set_alt(InstructionsList[-4])

#THE FITNESS FUNCTION
fit_ley <- function(s) {
             require("RecordLinkage")
             return(levenshteinSim(s,"baobab"))
           }

#THE HAPLOTYPE
MyHAP         <- Haplotype_class$new(HAP=sample(InstructionsList,6, replace=T),
                                         mutated=F
                                    )
#THE ENGINE
MyEngine      <- Engine_Class$new(population=list(),
                                  mutation.rate=0.1,
                                  fitness.function=fit_ley,
                                  step = 1L
                                  )

#THE ENGINE'S POPULATION
MyEngine$population <- lapply(rep(NA,200), function(x) {Individual_Class$new(
                                                          genome=MyHAP,
                                                          age=0L,min.fit=Inf,
                                                          fit.sum=0L,
                                                          local.engine=MyEngine
                                              )
                                           }
                             )

