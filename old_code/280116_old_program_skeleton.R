################################################################################
## REFERENCECLASSES DECLARATION SECTION
# The referenceclass used for wrapping Dummy instructions
Engine_Class <- setRefClass(Class='Engine',
                            fields=list(population    = 'list',
                                        step    = 'integer',
                                        mutation.rate = 'numeric',
                                        fitness.function = 'function'
                                       )
                           )

Haplotype_class <- setRefClass(Class='HAP',
                               fields=list(HAP    = 'list',
                               mutated = 'logical'
                              )
)

Individual_Class <- setRefClass(Class='Individual',
                                fields=list(genome       = 'HAP',
                                            phenotype       = 'list',
                                            fitness      = 'numeric',
                                            local.engine = 'Engine',
                                            min.fit      = 'numeric',
                                            fit.sum      = 'numeric',
                                            age          = 'numeric'
                                           )
)

Instruction_Class <- setRefClass(Class='Instruction',
                                 fields=list(instruction       = 'function',
                                             alternatives.list = 'list'
                                            )
)

Instruction_Class$methods(
  get_instr = function() {
               return(instruction)
              },
  set_instr = function(new.inst) {
               instruction <<- new.inst
              },
  run       = function(ind.run) {
               pos = ind.run$age
               ind.run$phenotype[[pos]] <- instruction(ind.run$phenotype[[pos]])
              },
  set_alt    = function(new.list) {
                alternatives.list <<- unlist(new.list)
               },
  show       = function () {
                              if (is.null(cl <- tryCatch(class(.self),
                                                         error = function(e)
                                                                  NULL
                                                        )
                                         )
                                 ) {
                                  cat("Prototypical reference class object\n")
                              }
                              else {
                                  cat("Reference class object of class ",
                                      classLabel(cl), 
                                      "\n", sep = "")
                                  fields <- names(.refClassDef@fieldClasses)
                                  for (fi in fields) {
                                      cat("Field \"", fi, "\":\n", sep = "")
								      if (typeof(field(fi)) == 'list') {
										print(names(field(fi)))
                                      } else {
                                          methods::show(field(fi))
                                        }
                                  }
                              }
                          },

  get_alt    = function() {
                out = sample(alternatives.list,1)
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
                      mutated <<- TRUE
                   },
  run            = function(individual) {
                    for (i in seq_along(HAP)) {
                     val=HAP[[i]]$run(individual)
                    }
                   }
)

Individual_Class$methods(
  get_ploidy     = function() {
                    return(length(.self$genome))
                   },
  get_genome     = function() {
                    return(.self$genome)
                   },
  set_genome     = function(HAP) {
                    .self$genome <- list(HAP)
                   },
  produce_gamete = function(.self,mutation.rate) {
                    offspring.genome <- NA
                    if(.self$get_ploidy() == 1){
                     genome.len <- length(.self$genome$get_HAP()) 
                     sid = runif(1)
                     mut.instr.num <- min(mutation.rate/sid,genome.len)
                     if (mut.instr.num == 0) {
                      offspring.genome <- .self$genome
                     } else {
                      offspring.genome <- .self$genome$copy(shallow=T)
                      mut.instr.idx = sample(length(offspring.genome$get_HAP()),
                                     mut.instr.num
                                    )
                      offspring.genome$mutate(mut.instr.idx)
                     }
                    } else {
                     print("produce gamete not implemented for ploidy > 1")
                    }
                    return(offspring.genome)
                   },
  get_results    = function(value=length(phenotype)) {
                    return(phenotype[[value]])
                   },
  reproduce      = function(){
                    gamete = produce_gamete(local.engine$mutation.rate)
                    if (length(gamete) == 1){
                     return(Individual_Class$new(genome=gamete,
                                                 min.fit=Inf,
                                                 fit.sum=0L,
                                                 age=0L,
                                                 local.engine = local.engine
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
  get_avg.fit    = function(step.n = .self$age) {
                              return(.self$fit.sum / .self$age)
                              },
  run            = function(.self) {
                    .self$age <- .self$age + 1L
                    position <- .self$age
                    .self$phenotype <- append(.self$phenotype,'')
                    .self$genome$run(.self)
                    .self$fitness <- append(.self$fitness,0)
                    .self$fitness[position] <- .self$local.engine$fitness.function(
                                                 .self$phenotype[[position]]
                                               )
                    if ( .self$fitness[position] < .self$min.fit){
                      .self$min.fit <- .self$fitness[position]
                    }
                    .self$fit.sum <- .self$fit.sum + .self$fitness[position]
                    return()
                   }
)

#The referenceclass used of the engine

Engine_Class$methods(
  add_organism = function(.self, ind) {
                  ind$attach_engine(.self)
                  population <<- append(population,ind)
                 },
  run          = function(.self, stepsA) {
                  for (i in seq_len(stepsA)){
                   sapply(.self$population,
                          function(x) {x$run()}
                         )
                  if (.self$step %% 2 == 0){
#effective.population <- max(1,round(length(.self$population)/10))
#reproducers <- sample(.self$population,effective.population)
#newborns    <- lapply(reproducers,
                  .self$prune()
                     newborns    <- lapply(.self$population,
                                          function(x){x$reproduce()})
                     .self$population <- append(newborns,.self$population)
                   }
                  .self$step <-.self$ step + 1L
                  pop.fit.lis <- vapply(.self$population, function(x) {
                                                           x$get_avg.fit()
                                                          }
                                                        , FUN.VALUE=1)
                  pop.fit.lis[is.nan(pop.fit.lis)] <- NA
                  avg.pop.fit <- mean(unlist(pop.fit.lis), na.rm=T)
                  cat(.self$step," | ",
                      length(.self$population),' |',
                      avg.pop.fit,' | ',
                      .self$get_ref(), '\n'
                     )
                  }
                 },
  restart      = function(){
                  step <<- 0
                 },
  set_fitness.function = function(fit.fun.out){
                           fitness.function <<- fit.fun.out
                         },
  prune        = function(){
                  k.fit.ord <- vapply(population, function(x){
                                             #x$fitness[x$age]
                                             x$get_min.fit()
                                            },
                                      FUN.VALUE=1
                                     )
                                     #print(k.fit.ord)
                  population <<- population[order(k.fit.ord, decreasing = T)]
                  population <<- population[1:(length(population)/2)]
                 },
  get_ref = function(){
                        table(unlist(lapply(population, function(x){x$phenotype[1]})))
  }
)

# Dummy list of instructions
DummyList <- list(
#DummyInstruction0 = function(x) {return(x)},
DummyInstruction1 = function(x) {return(paste(x, 'b', sep=''))},
DummyInstruction2 = function(x) {return(paste(x, 'a', sep=''))},
DummyInstruction3 = function(x) {return(paste(x, 'o', sep=''))}
)
# The list containing wrappers to the fundamental instructions
InstructionsList <- lapply(DummyList, function(x) {Instruction_Class$new(instruction=x)})
#names <- paste ('Fun',1:3, sep='')
#invisible(mapply(function(x,y) { x$name <- y}, InstructionsList,names))
names(InstructionsList) = paste('Instruction', seq_along(DummyList), sep='')
InstructionsList[[1]]$set_alt(InstructionsList[-1])
InstructionsList[[2]]$set_alt(InstructionsList[-2])
InstructionsList[[3]]$set_alt(InstructionsList[-3])
#InstructionsList[[4]]$set_alt(InstructionsList[-4])
#The fitness function
fit_ley <- function(s) {
             require('RecordLinkage')
             return(levenshteinSim(s,'baobab'))
           }

MyHAP         <- Haplotype_class$new(HAP=sample(InstructionsList,6, replace=T), mutated=F)
MyEngine      <- Engine_Class$new(population=list(),
                                   mutation.rate=0.3,
                                  fitness.function=fit_ley,
                                  step = 1L
                                  )

MyEngine$population <- lapply(rep(NA,200), function(x) {Individual_Class$new(genome=MyHAP, age=0L,min.fit=Inf,fit.sum=0L,local.engine=MyEngine)})










#MyIndividual  <- Individual_Class$new(genome=MyHAP,
#                                      age    =0L,
#                                      min.fit=0L,
#                                      fit.sum=0L
#                 )
#w,,ith(MyIndividual, list(phenotype,fitness,min.fit,fit.sum))
# table(unlist(lapply(MyEngine$population, function(x){x$phenotype[1]})))
#xsel | awk -F "|" '{print $4}' | column -t -s ' ' | vim -

# add method in simiand to create n individuals
#store parameters in vengine
#add method in engine to return pop fit
