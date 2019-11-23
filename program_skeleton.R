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
#                                             args    = "list"
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
#  get_arg  = function() {
#               return(arg)
#             },
#  set_arg  = function(...) {
#               arg <<- unlist(...)
#             },
  mutate   = function(sites) {
                HAP[sites] <<- lapply(HAP[sites], function(x) {
                                                               x <- x$get_alt()           
                                                             }
                                          )
             mutated <<- TRUE
             },
  run      = function(individual) {
               lapply(HAP, function(x){x$run(individual)})
#               val = mapply(function(x,y) {
#                              ili = y[1:length(formals(x))]
#                              do.call(x,ili)
#                            }
#                            HAP,
#                            arg
#                           )
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
                     genome <<- list(HAP)
                   },
  produce_gamete = function(mutation.rate) {
                     offspring.genome <- NA
                     if (get_ploidy() == 1){
                       genome.len <- length(genome$get_HAP())
                       sid = runif(1)
                       mut.instr.num <- min(mutation.rate/sid,genome.len)
                       if (mut.instr.num == 0) {
                         offspring.genome <- genome
                       } else {
                           offspring.genome <- genome$copy(shallow=T)
                           #offspring.genome <- genome$copy()
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
                     return(fit.sum / age)
                     },
  run            = function() {
                     age <<- age + 1L
                     position <- age
                     phenotype <<- append(phenotype,"")
                     genome$run(.self)
                     fitness <<- append(fitness,0)
                     fitness[position] <<- local.engine$fitness.function(
                                                  phenotype[[position]]
                                                )
                     if (fitness[position] < min.fit){
                       min.fit <<- fitness[position]
                     }
                     fit.sum <<- fit.sum + fitness[position]
                   }
)

Engine_Class$methods(
  add_organism = function(ind) {
                   ind$attach_engine(.self)
                   population <<- append(population,ind)
                 },
  run          = function(stepsA) {
                   for (i in seq_len(stepsA)){
                     sapply(population, function(x) {x$run()})
                     pop.size <- length(population)
                     if (step %% 2 == 0){
                     pop.min.fit <- vapply(population, function(x){
                                                         x$fitness[1]
                                                       },
                                           FUN.VALUE=1
                                          )
                     k.fit.ord <- rank(pop.min.fit, ties.method="first")
                     are.fit <- k.fit.ord > 1/2 * pop.size
                     will.be.replaced <- k.fit.ord <= 1/2 * pop.size
                     population[will.be.replaced] <<- lapply(
                                                             population[are.fit],
                                                             function(x){
                                                               x$reproduce()
                                                             }
                                                             )
                     }
                     step <<-step + 1L
                     print(step)
                     #print(cbind(fit_ley(names(.self$get_ref())),
                     #      .self$get_ref()
                     #     ))
                     #pop.fit.lis <- vapply(population, function(x) {
                     #                                          x$get_avg.fit()
                     #                                        }
                     #                                      , FUN.VALUE=1
                     #                     )
                     #pop.fit.lis[is.nan(pop.fit.lis)] <- NA
                     #avg.pop.fit <- mean(unlist(pop.fit.lis), na.rm=T)
                     #cat(step," | ", length(population)," |",
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
