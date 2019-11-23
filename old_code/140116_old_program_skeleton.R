################################################################################
Engine_Class <- setRefClass(Class  = "Engine",
                            fields = list(population       = "list",
                                          step             = "integer",
                                          mutation.rate    = "numeric",
                                          fitness.function = "function"
                                         )
                           )

Haplotype_class <- setRefClass(Class   = "HAP",
                               fields  = list(HAP     = "list",
                                              mutated = "logical"
                              )
)

Individual_Class <- setRefClass(Class  = "Individual",
                                fields = list(genome       = "HAP",
                                              phenotype    = "list",
                                              fitness      = "numeric",
                                              local.engine = "Engine",
                                              min.fit      = "numeric",
                                              fit.sum      = "numeric",
                                              age          = "integer"
                                             )
)

Instruction_Class <- setRefClass(Class  = "Instruction",
                                 fields = list(instruction       = "function",
                                               alternatives.list = "list"
                                              )
)


Instruction_Class$methods(
  get_inst = function() {
               return(instruction)
             },
  set_inst = function(new.inst) {
               instruction <<- new.inst
             },
  run      = function(ind.run) {
               pos = ind.run$age
               ind.run$phenotype[[pos]] <- instruction(ind.run$phenotype[[pos]])
             },
  set_alt   = function(new.list) {
                alternatives.list <<- unlist(new.list)
              },
  show      = function () {
                if (is.null(cl <- tryCatch(class(.self),
                                           error = function(e) {NULL}
                                           )
                           )
                   ) {
                  cat("Prototypical reference class object\n")
                } else {
                        cat("Reference class object of class ",
                            classLabel(cl), "\n", sep = "")
                        fields <- names(.refClassDef@fieldClasses)
                        for (fi in fields) {
                          cat("Field \"", fi, "\":\n", sep = "")
		        	       if (typeof(field(fi)) == "list") {
		        	         print(names(field(fi)))
                          } else {
                              methods::show(field(fi))
                          }
                        }
                  }
              },
  get_alt   = function() {
                out = sample(alternatives.list, 1)
                return(out[[1]])
              }
)

Haplotype_class$methods(
  get_HAP  = function() {
               return(HAP)
             },
  set_HAP  = function(...) {
               HAP <<- unlist(...)
             },
  mutate   = function(sites) {
               HAP <<- lapply(HAP, function(x){x <- x$get_alt()})
               mutated <<- TRUE
             },
  run      = function(individual) {
               val = lapply(HAP, function(x){x$run(individual)})
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
                     if (.self$get_ploidy() == 1){
                       genome.len <- length(.self$genome$get_HAP())
                       sid = runif(1)
                       mut.instr.num <- min(mutation.rate/sid,genome.len)
                       if (mut.instr.num == 0) {
                         offspring.genome <- .self$genome
                       } else {
                           offspring.genome <- .self$genome$copy(shallow=T)
                           offspring.gen.len = length(offspring.genome$get_HAP()
                                                     )
                           mut.instr.idx = sample(offspring.gen.len,
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
                                                  )
                            )
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
                     .self$phenotype <- append(.self$phenotype,"")
                     .self$genome$run(.self)
                     .self$fitness <- append(.self$fitness,0)
                     .self$fitness[position] <- .self$local.engine$fitness.function(
                                                  .self$phenotype[[position]]
                                                )
                     if (.self$fitness[position] < .self$min.fit){
                       .self$min.fit <- .self$fitness[position]
                     }
                     .self$fit.sum <- .self$fit.sum + .self$fitness[position]
                   }
)

Engine_Class$methods(
  add_organism = function(.self, ind) {
                   ind$attach_engine(.self)
                   population <<- append(population,ind)
                 },
  run          = function(.self, stepsA) {
                   for (i in seq_len(stepsA)){
                     sapply(.self$population, function(x) {x$run()})
                     if (.self$step %% 2 == 0){
                       #effective.population <- max(1,round(length(.self$population)/10))
                       #reproducers <- sample(.self$population,effective.population)
                       #newborns    <- lapply(reproducers,
                       #.self$prune()
                  k.fit.ord <- vapply(.self$population, function(x){ x$get_min.fit()},
                                      FUN.VALUE=1
                                     )
                  .self$population <- .self$population[order(k.fit.ord, decreasing = T)]
                       reproducers <- .self$population[1:(length(.self$population)/2)] 
                       newborns    <- lapply(reproducers, function(x){
                                                                 x$reproduce()
                                                               }
                                            )
                  .self$population[!(.self$population %in% reproducers)]<- newborns 
                     }
                     .self$step <-.self$ step + 1L
                     #pop.fit.lis <- vapply(.self$population, function(x) {
                     #                                          x$get_avg.fit()
                     #                                        }
                     #                                      , FUN.VALUE=1
                     #                     )
                     #pop.fit.lis[is.nan(pop.fit.lis)] <- NA
                     #avg.pop.fit <- mean(unlist(pop.fit.lis), na.rm=T)
                     #cat(.self$step," | ", length(.self$population)," |",
                     #    avg.pop.fit," | ", .self$get_ref(), "\n"
                     #   )
                   }
                 },
  restart      = function(){
                   step <<- 0
                 },
  set_fitness.function = function(fit.fun.out){
                           fitness.function <<- fit.fun.out
                         },
  prune        = function(){
                  k.fit.ord <- vapply(population, function(x){ x$get_min.fit()},
                                      FUN.VALUE=1
                                     )
                  population <<- population[order(k.fit.ord, decreasing = T)]
                  population <<- population[1:(length(population)/2)]
                 },
  get_ref = function(){
              table(unlist(lapply(population, function(x){x$phenotype[1]})))
            }
)
