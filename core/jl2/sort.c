#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include "name.h"
#include "errmem.h"
#include "types.h"

#define T unsigned int
#define SORT_SUFFIX _ui
#include "sort_imp.c"
#undef SORT_SUFFIX
#undef T
#if 1
#define T unsigned long
#define SORT_SUFFIX _ul
#include "sort_imp.c"
#undef SORT_SUFFIX
#undef T
#endif
