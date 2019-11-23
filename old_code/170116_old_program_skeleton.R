##DECLARATION SECTION
#The number of fundamental instructions
FundamentalInstructionsNumber <- 100
#The environment where the instruction are stored
InstructionEnvironment <- new.env()
#The environment where the instruction Classs are stored
InstructionClassEnvironment <- new.env()
#The length of the main function
MainFunctionLength=10

##DEFINITION OF THE REFERENCECLASS USED FOR STORING INSTRUCTIONS
InstructionClass <- setRefClass(Class='InstructionClass', fields=list(instruction='function'))
HAP_Class <- setRefClass(Class='HAP_Class', fields=list(HAP='list')) #ADDED 1701

##POPULATE THE InstructionEnvironment WITH DUMMY FUNCTIONS
for (i in seq(FundamentalInstructionsNumber)) { assign(paste('instruction', i, sep=''), function(x) {return(x)}, envir=InstructionEnvironment) }


##POPULATE THE InstructionClassEnvironment WITH InstructionClass objects
#(one per instruction)
instructionList <- as.list(InstructionEnvironment)
X <- lapply(instructionList, function(x) InstructionClass$new(instruction=x))
names(X) <- paste('InstructionClass', seq(100), sep='')
InstructionClassEnvironment <<- list2env(X)
rm(instructionList, X, i)
#REMOVED 1701 for (i in seq(along=InstructionEnvironment)){
#REMOVED 1701     assign(paste('InstructionClass', i, sep=''),
#REMOVED 1701         InstructionClass$new(instruction=instructionList[[i]]),
#REMOVED 1701         envir=InstructionClassEnvironment
#REMOVED 1701         )
#REMOVED 1701     }


##CREATE THE MAIN FUNCTION
#namely a list of random pointers
#sampled from InstructionClassEnvironment
HAP_1 <- HAP_Class$new(HAP=list(InstructionClassEnvironment$InstructionClass1, InstructionClassEnvironment$InstructionClass2)) #ADDED 1801
#REMOVED 1801 MainFunction <- list(sample(as.list(InstructionClassEnvironment), MainFunctionLength))
