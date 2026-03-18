#ifndef MAIN_H
#define MAIN_H

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <unordered_map>
#include <filesystem>
#include <vector>
#include <istream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include <getopt.h>

#include "log.h"
#include "uid-generator.h"
#include "bed.h"
#include "global.h"
#include "struct.h"
#include "functions.h"

#include "gfa-lines.h"
#include "gfa.h"
#include "sak.h"

#include "stream-obj.h"

#include "input-agp.h"
#include "input-filters.h"
#include "input-gfa.h"

#include "input.h"
#include "tools.h"

#endif /* MAIN_H */
