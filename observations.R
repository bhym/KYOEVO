#MANUALLY RUN EACH POPULATION
sapply(population, function(x) x$run())

#MANUALLY DO SELECTION
pop.min.fit <- vapply(population, function(x){ x$get_min.fit()}, FUN.VALUE=1)
k.fit.ord <- order(pop.min.fit, decreasing=T)
are.fit <- k.fit.ord > 1/2 * pop.size
will.be.replaced <- k.fit.ord <= 1/2 * pop.size
population[will.be.replaced] <- lapply(population[are.fit], function(x){x$reproduce()})
print(cbind(fit_ley(names(table(unlist(lapply(population,
function(x){x$phenotype[1]}))))),table(unlist(lapply(population,
function(x){x$phenotype[1]}))))
)


#with(MyIndividual, list(phenotype,fitness,min.fit,fit.sum))

##COUNT TABLE OF PHENOTYPES
#table(unlist(lapply(MyEngine$population, function(x){x$phenotype[x$age]})))

##TABLE OF PHENOTYPE,FITNESS AND MIN FITNESS. ROWS ARE THE INDIVIDUALS
#t(sapply(MyEngine$population, function(x){cbind(x$phenotype[1],x$fitness[1],x$min.fit[1])}))

#LAST THING I DID
#add method in simiand to create n individuals
#store parameters in vengine
#Implement function arguments.
#In this way, if the user-defined R-function requires arguments,
#it can get them by calling organism$get_HAP()$get_arg(pos).

