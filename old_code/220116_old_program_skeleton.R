## REFERENCECLASSES DECLARATION SECTION
# The referenceclass used for wrapping Dummy instructions
Instruction_Class <- setRefClass(Class='Instruction', fields=list(instr   = 'function',
                                                                  altList = 'list'))
Instruction_Class$methods(get_instr = function() {
                                       return(instr)
                                      },
                          set_instr = function(inst) {
                                       instr <<- inst
                                      },
                          run       = function(ind) {
                                       ind$output[[ind$iter]] <- instr(ind$output[[ind$iter]])
                                      },
                         get_alt    = function() {
                                       out = sample(altList,1)
                                       return(out)
                                      }
                         )

# The referenceclass used for Haplotypes (i.e. functions)
Haplotype_class <- setRefClass(Class='HAP', fields=list(HAP    = 'list',
                                                        mutatd = 'logical')
                        )
Haplotype_class$methods(get_HAP  = function() {
                                    return(HAP)
                                   },
                  set_HAP        = function(...) {
                                    HAP <<- unlist(...)
                                   },
                  mutate         = function(mutate, pos) {
                                     if (mutate==T) {
                                      HAP[[pos]] <<- HAP[[pos]]$get_alt()
                                      mutatd <<- TRUE
                                     } else {
                                      mutatd <<- FALSE
                                     }
                                   },
                  run            = function(ind) {
                                    for (i in seq(along=HAP)) {
                                     val=HAP[[i]]$run(ind)
                                    }
                                   }
                )

#The referenceclass used for individuals
Individual_Class <- setRefClass(Class='Individual', fields=list(genome = 'list',
                                                                iter   = 'integer',
                                                                output = 'list')
                                                                )
Individual_Class$methods(get_ploidy     = function() {
                                           return(length(genome))
                                          },
                         get_genome     = function() {
                                           return(genome)
                                          },
                         set_genome     = function(HAP) {
                                           genome <<- list(HAP)
                                          },
                         produce_gamete = function() {
                                           return(genome[[sample(length(genome),1)]]$HAP)
                                          },
                         get_results    = function(value=length(output)) {
                                           return(output[[value]])
                                          },
                         #mutated        = function() {
                         #                  mutcheck = sapply(genome,
                         #                                    function(x) {x$mutatd}
                         #                                   )
                         #                  check    = any(mutcheck)
                         #                  return(check)
                         #                 },
                         run            = function(.self) {
                                           .self$iter <- .self$iter +1L
                                           .self$output <- append(.self$output,'')
                                           sapply(.self$genome,
                                                  function(x) {x$run(.self)}
                                                 )
                                           return()
                                          }
                        )

#The referenceclass used of the engine
Engine_Class <- setRefClass(Class='Engine', fields=list(inds      = 'list',
                                                        iteration = 'integer',
                                                        mutrate   = 'numeric'
                                                        )
                            )

Engine_Class$methods(add_organism = function(ind) {
                                     inds <<- append(inds,ind)
                                    },
                     run          = function(steps) {
                                     for (i in seq_len(steps)){
                                      #sapply(inds,
                                      #       function(x) {x$mutated(mutrate)}
                                      #      )
                                      sapply(inds,
                                             function(x) {x$run()}
                                            )
                                     }
                                    },
                     restart      = function(){
                                     iteration <<- 0
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
MyHAP         <- Haplotype_class$new(HAP=InstructionsList[1:2], mutatd=F)
MyIndividual  <- Individual_Class$new(genome=list(MyHAP), iter=0L)
MyEngine      <- Engine_Class$new(inds=list(MyIndividual), iteration=0L, mutrate=0.75e-4)
