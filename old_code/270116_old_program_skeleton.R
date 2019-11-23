################################################################################
## REFERENCECLASSES DECLARATION SECTION
# The referenceclass used for wrapping Dummy instructions
Engine_Class <- setRefClass(Class='Engine', fields=list(inds      = 'list',
                                                        generation = 'integer',
                                                        mutrate   = 'numeric',
                                                        fitness.function = 'function'
                                                        )
                            )

Individual_Class <- setRefClass(Class='Individual',
                                fields=list(genome = 'HAP',
                                            iter   = 'integer',
                                            output = 'list',
                                            fitness= 'numeric',
                                            local.engine = 'Engine',
                                            min.fit = 'numeric',
                                            fit.sum = 'numeric'
                                           )
)

Haplotype_class <- setRefClass(Class='HAP', fields=list(HAP    = 'list',
                                                        mutatd = 'logical')
)
Instruction_Class <- setRefClass(Class='Instruction',
                                 fields=list(instr   = 'function',
                                             altList = 'list'
                                            )
)

Instruction_Class$methods(
  get_instr = function() {
               return(instr)
              },
  set_instr = function(inst) {
               instr <<- inst
              },
  run       = function(ind) {
               ind$output[[ind$iter]] <- instr(ind$output[[ind$iter]])
              },
  set_alt    = function(altlist) {
                altList <<- unlist(altlist)
               },
  show       = function () {
                              if (is.null(cl <- tryCatch(class(.self), error = function(e) NULL))) {
                                  cat("Prototypical reference class object\n")
                              }
                              else {
                                  cat("Reference class object of class ", classLabel(cl), 
                                      "\n", sep = "")
                                  fields <- names(.refClassDef@fieldClasses)
                                  for (fi in fields) {
                                      cat("Field \"", fi, "\":\n", sep = "")
								      if (typeof(field(fi)) == 'list') {
										print(names(field(fi)))
                                      } else {
                                          methods::show(field(fi))
                                        }
									  #print(class(fields))
								      #print(typeof(field(fi)[1]))
                                  }
                              }
                          },

  get_alt    = function() {
               out = sample(altList,1)
               return(out[[1]])
              }
)

Haplotype_class$methods(
  get_HAP  = function() {
                    return(HAP)
                   },
  set_HAP        = function(...) {
                    HAP <<- unlist(...)
                   },
  mutate         = function(sites) {
                      for (i in sites){
                       HAP[[i]] <<- HAP[[i]]$get_alt()
                      }
                      mutatd <<- TRUE
                   },
  run            = function(ind) {
                    for (i in seq_along(HAP)) {
                     val=HAP[[i]]$run(ind)
                    }
                   }
)

Individual_Class$methods(
  get_ploidy     = function(.self) {
                    return(length(.self$genome))
                   },
  get_genome     = function() {
                    return(.self$genome)
                   },
  set_genome     = function(HAP) {
                    .self$genome <- list(HAP)
                   },
  produce_gamete = function(.self,rate) {
                    offspring_genome <- NA
                    if(.self$get_ploidy() == 1){
                     nb <- length(.self$genome$get_HAP()) * rate
                     if (nb == 0) {
                      offspring_genome <- .self$genome
                     } else {
                      offspring_genome <- .self$genome$copy(shallow=T)
                      ogl              <- length(offspring_genome$HAP)
                      pippo = sample(ogl,nb)
                      offspring_genome$mutate(pippo)
                     }
                    } else {
                     print("produce gamete not implemented for ploidy > 1")
                    }
                    return(offspring_genome)
                   },
  get_results    = function(value=length(output)) {
                    return(output[[value]])
                   },
  reproduce      = function(){
                    gamete = produce_gamete(local.engine$mutrate)
                    if (length(gamete) == 1){
                     return(Individual_Class$new(genome=gamete,
                                                 iter=0L,
                                                 min.fit=0L,
                                                 fit.sum=0L,
                                                 local.engine=local.engine
                                                ))
                    } else {
                     print("Reproduction not implemented for ploidy >1")
                    }
                   },
  attach_engine  = function(engine.to.attach){
                              local.engine <<- engine.to.attach
                   },
  get_min.fit    = function() {
                              return(min.fit)
                              },
  get_avg.fit    = function(i) {
                              return(fitness[i] / fit.sum)
                              },
  run            = function(.self) {
                    .self$iter <- .self$iter +1L
                    .self$output <- append(.self$output,'')
                    .self$genome$run(.self)
                    .self$fitness <- append(.self$fitness,0)
                    .self$fitness[.self$iter] <- .self$local.engine$fitness.function(.self$output[[.self$iter]])
                    if ( .self$fitness[.self$iter] < .self$min.fit){
                      .self$min.fit <- .self$fitness[.self$iter]
                    }
                    .self$fit.sum <- .self$fit.sum + .self$fitness[.self$iter]
                    return()
                   }
)

#The referenceclass used of the engine

Engine_Class$methods(
  add_organism = function(.self, ind) {
                  ind$attach_engine(.self)
                  inds <<- append(inds,ind)
                 },
  run          = function(.self, steps) {
                  .self$prune()
                  for (i in seq_len(steps)){
                   sapply(.self$inds,
                          function(x) {x$run()}
                         )
                  if (.self$generation %% 11 == 0){
                      reproducers <- sample(.self$inds,max(1,round(length(.self$inds)/10)))
                     newborns    <- lapply(reproducers,
                                          function(x){x$reproduce()})
                     .self$inds <- append(newborns,.self$inds)
                   }
                  .self$generation <-.self$ generation + 1L
                  print(.self$generation)
                  }
                 },
  restart      = function(){
                  generation <<- 0
                 },
  set_fitness.function = function(fit.fun){
                           fitness.function <<- fit.fun
                         },
  prune        = function(){
                  k.fit.ord <- vapply(inds, function(x){x$get_min.fit()}, FUN.VALUE=1)
                  inds <<- inds[order(k.fit.ord)]
                  inds <<- inds[1:(length(inds)/2)]
                 }
)

# Dummy list of instructions
DummyList <- list(
DummyInstruction0 = function(x) {return(x)},
DummyInstruction1 = function(x) {return(paste(x, 'b', sep=''))},
DummyInstruction2 = function(x) {return(paste(x, 'a', sep=''))},
DummyInstruction3 = function(x) {return(paste(x, 'o', sep=''))},
DummyInstruction4 = function(x) {return(paste(x, ' ', sep=''))}
)
# The list containing wrappers to the fundamental instructions
InstructionsList <- lapply(DummyList, function(x) {Instruction_Class$new(instr=x)})
#names <- paste ('Fun',1:3, sep='')
#invisible(mapply(function(x,y) { x$name <- y}, InstructionsList,names))
names(InstructionsList) = paste('Instruction', seq_along(DummyList), sep='')
InstructionsList[[1]]$set_alt(sample(InstructionsList,2))
InstructionsList[[2]]$set_alt(sample(InstructionsList,2))
InstructionsList[[3]]$set_alt(sample(InstructionsList,2))
InstructionsList[[4]]$set_alt(sample(InstructionsList,2))
InstructionsList[[5]]$set_alt(sample(InstructionsList,2))
#The fitness function
fit_ley <- function(s) {
             require('RecordLinkage')
             return(levenshteinSim(s,'baobab'))
           }

MyHAP         <- Haplotype_class$new(HAP=sample(InstructionsList,6, replace=T), mutatd=F)
MyIndividual  <- Individual_Class$new(genome=MyHAP,
                                      iter=0L,
                                      min.fit=0L,
                                      fit.sum=0L
                 )
MyEngine      <- Engine_Class$new(inds=list(MyIndividual),
                                   mutrate=1,
                                  fitness.function=fit_ley,
                                  generation = 1L
                                  )


#with(MyIndividual, list(iter,output,fitness,min.fit,fit.sum))
