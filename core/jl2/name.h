#ifndef NAME_H
#define NAME_H

/* establishes some macros to establish
   * the FORTRAN naming convention
     default      gs_setup, etc.
     -DUPCASE     GS_SETUP, etc.
     -DUNDERSCORE gs_setup_, etc.
   * a prefix for all external (non-FORTRAN) function names
     for example, -DPREFIX=jl_   transforms fail -> jl_fail
   * a prefix for all external FORTRAN function names     
     for example, -DFPREFIX=jlf_ transforms gs_setup_ -> jlf_gs_setup_
*/

/* the following macro functions like a##b,
   but will expand a and/or b if they are themselves macros */
#define TOKEN_PASTE_(a,b) a##b
#define TOKEN_PASTE(a,b) TOKEN_PASTE_(a,b)

#if defined(FPREFIX)
#  if defined(UPCASE)
#    define FORTRAN_NAME(low,up) TOKEN_PASTE(FPREFIX,up)
#  elif defined(UNDERSCORE)
#    define FORTRAN_NAME(low,up) TOKEN_PASTE(TOKEN_PASTE(FPREFIX,low),_)
#  else
#    define FORTRAN_NAME(low,up) TOKEN_PASTE(FPREFIX,low)
#  endif
#else
#  if defined(UPCASE)
#    define FORTRAN_NAME(low,up) up
#  elif defined(UNDERSCORE)
#    define FORTRAN_NAME(low,up) TOKEN_PASTE(low,_)
#  else
#    define FORTRAN_NAME(low,up) low
#  endif
#endif

#endif

