# KYOEVO

This repository contains the code I wrote together with Prof. Chris Vincenot during my visit as PhD student to the biosphere lab of the Dept. of Social Informatics within the Kyoto University.  The code is working, but in alpha stage: at the moment it correctly initialises an environment, and populates it with artificial programs able to execute simple operations. A fitness function evaluates the operation outputs, and selects the organisms with the best one.

## To do
* #### Add function arguments
  * Implement arithmetics game (i.e. using operators +,-,/,*, and numbers from 0 to 9, find calculation closest to 137)

* #### Add individual variables
  * Add "()" to arithmetics game. Done by having a local variable which would be a stack. "(" -> push, ")"->pull, and expressions are right-handed. e.g. (4 5 +)(3 1 -)*
  * Think about IBM capabilities

* #### Add environmental variables
  * Starting from 1 organism, try to evolve two mutualistic species that "collaborate" to reach a goal. This could be to generate a string (or maybe even to play the arithmetics game). Individuals have a finite code length, they take strings from the environment, and drop the resulting output string back in the environment. They are selected based on progress (difference between input and output string distance from the target). Open question: how to manage the environment string pool?  

* #### Add B-functions (as opposed to R-functions).
  * These are sets of instructions.
