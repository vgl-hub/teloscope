#include <stdlib.h>
#include <string>

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "stream-obj.h"
#include "fastx.h"

#include "gfa-lines.h"
#include "gfa.h"

#include "teloscope.h"
#include "input.h"

void Input::load(UserInput userInput) {
    
    this->userInput = userInput;
    
}

void Input::read(InSequences &inSequence) {

    loadSequences(userInput, &inSequence, 'f', 0); // load from FASTA/FASTQ to templated object

    jobWait(threadPool);


}