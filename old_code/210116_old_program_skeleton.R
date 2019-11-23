## REFERENCECLASSES DECLARATION SECTION
# The referenceclass used for wrapping Dummy instructions
Instruction_Class <- setRefClass(Class='Instruction', fields=list(instr='function'))
Instruction_Class$methods(get_instr = function() {
                                       return(instr)
                                      },
                          set_instr = function(inst) {
                                       instr <<- inst
                                      },
                          run       = function(j) {
                                       return(instr(j))
                                      }
                         )

# The referenceclass used for Haplotypes (i.e. functions)
HAP_Class <- setRefClass(Class='HAP', fields=list(HAP    = 'list',
                                                  mutatd = 'logical')
                        )
HAP_Class$methods(get_HAP = function() {
                             return(HAP)
                            },
                  set_HAP = function(...){
                             HAP <<- unlist(...)
                            },
                  mutate  = function(prob, pos, InstList) {
                              if (prob==T) {
                               new_inst <- sample(length(InstList),1)
                               HAP[[pos]] <<- InstList[[new_inst]]
                               mutatd <<- TRUE
                              } else {
                               mutatd <<- FALSE
                              }
                            },
                  run     = function(init) {
                             val=init
                             for (i in seq(along=HAP)) {
                              val=HAP[[i]]$run(val)
                             }
                             return(val)
                            }
                )

#The referenceclass used for individuals
Individual_Class <- setRefClass(Class='Individual', fields=list(genome='list',
                                                                iter  ='integer',
                                                                output='list')
                                                                )
Individual_Class$methods(get_ploidy  = function() {
                                        return(length(genome))
                                       },
                         get_genome  = function() {
                                        return(genome)
                                       },
                         set_genome  = function(HAP) {
                                        genome <<- list(HAP)
                                       },
                         get_results = function(value=length(output)) {
                                        return(output[[value]])
                                       },
                         mutated     = function() {
                                        mutcheck = sapply(genome,
                                                          function(x) {x$mutatd}
                                                         )
                                        check    = any(mutcheck)
                                        return(check)
                                       },
                         run         = function(sv) {
                                        sv      = output[[length(output)]]
                                        cur_res = sapply(genome,
                                                         function(x) {x$run(sv)}
                                                        )
                                        output <<- append(output,cur_res)
                                       }
                        )

# Dummy list of instructions
DummyList <- list(
DummyInstruction1 = function(x) {return(paste(x, 1, sep=''))},
DummyInstruction2 = function(x) {return(paste(x, 2, sep=''))},
DummyInstruction3 = function(x) {return(paste(x, 3, sep=''))}
)
# The list containing wrappers to the fundamental instructions
InstructionsList <- lapply(DummyList, function(x) Instruction_Class$new(instr=x))
names(InstructionsList) = paste('Instruction', seq(along=DummyList), sep='')


#Three objects generated for testing the code
MyInstruction <- Instruction_Class$new(instr=function(x) {return('BOOM')})
MyHAP         <- HAP_Class$new(HAP=InstructionsList[1:2], mutatd=F)
MyIndividual  <- Individual_Class$new(genome=list(MyHAP), iter=0L, output=list(1))
