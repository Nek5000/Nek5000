/*
 *  $Id: mpi.h,v 1.45 1995/05/09 17:44:52 gropp Exp $
 *
 *  (C) 1993 by Argonne National Laboratory and Mississipi State University.
 *      All rights reserved.  See COPYRIGHT in top-level directory.
 */

/* user include file for MPI programs */

#ifndef _MPI_INCLUDE
#define _MPI_INCLUDE

/* Keep C++ compilers from getting confused */
#if defined(__cplusplus)
extern "C" {
#endif


/* Results of the compare operations */
/* These should stay ordered */
#define MPI_IDENT     0
#define MPI_CONGRUENT 1
#define MPI_SIMILAR   2
#define MPI_UNEQUAL   3

/* Data types */
/* MPI_LONG_LONG_INT is an OPTIONAL type, and may be set to NULL */
typedef struct MPIR_DATATYPE *MPI_Datatype;

/* 
   To allow compile-time use of MPI 'constants', they are declared in
   static storage.  This is incomplete at present. 

   Note that this will work for C but not for Fortran; in the Fortran
   case we will probably need predefined indices.  A C version that
   uses indices might have lower latency because of few cache misses
   while searching out the features of the basic datatypes.
 */
extern struct MPIR_DATATYPE MPIR_I_CHAR, MPIR_I_SHORT, MPIR_I_INT, MPIR_I_LONG,
                            MPIR_I_UCHAR, MPIR_I_USHORT, MPIR_I_UINT, 
                            MPIR_I_ULONG, MPIR_I_FLOAT, MPIR_I_DOUBLE, 
                            MPIR_I_LONG_DOUBLE, MPIR_I_LONG_DOUBLE_INT,
                            MPIR_I_BYTE, MPIR_I_PACKED, MPIR_I_UB, MPIR_I_LB,
                            MPIR_I_LONG_LONG_INT, MPIR_I_2INTEGER, 
                            MPIR_I_FLOAT_INT, MPIR_I_DOUBLE_INT, 
                            MPIR_I_LONG_INT, MPIR_I_SHORT_INT, MPIR_I_2INT,
                            MPIR_I_REAL, MPIR_I_DOUBLE_PRECISION, 
                            MPIR_I_COMPLEX, MPIR_I_DCOMPLEX, 
                            MPIR_I_LONG_DOUBLE, MPIR_I_LONG_LONG_INT, 
                            MPIR_I_LOGICAL;
#define MPI_CHAR (&MPIR_I_CHAR)
#define MPI_BYTE (&MPIR_I_BYTE)
#define MPI_SHORT (&MPIR_I_SHORT)
#define MPI_INT (&MPIR_I_INT)
#define MPI_LONG (&MPIR_I_LONG)
#define MPI_FLOAT (&MPIR_I_FLOAT)
#define MPI_DOUBLE (&MPIR_I_DOUBLE)

#define MPI_UNSIGNED_CHAR (&MPIR_I_UCHAR)
#define MPI_UNSIGNED_SHORT (&MPIR_I_USHORT)
#define MPI_UNSIGNED (&MPIR_I_UINT)
#define MPI_UNSIGNED_LONG (&MPIR_I_ULONG)

#define MPI_PACKED (&MPIR_I_PACKED)
#define MPI_UB (&MPIR_I_UB)
#define MPI_LB (&MPIR_I_LB)

#define MPI_FLOAT_INT (&MPIR_I_FLOAT_INT)
#define MPI_LONG_INT  (&MPIR_I_LONG_INT)
#define MPI_DOUBLE_INT (&MPIR_I_DOUBLE_INT)
#define MPI_SHORT_INT  (&MPIR_I_SHORT_INT)
#define MPI_2INT       (&MPIR_I_2INT)

extern MPI_Datatype MPI_LONG_DOUBLE, MPI_LONG_DOUBLE_INT, 
       MPI_LONG_LONG_INT, MPI_2INTEGER;
/* extern MPI_Datatype MPIR_complex_dte, MPIR_dcomplex_dte, 
   MPIR_logical_dte, ; */
extern MPI_Datatype MPIR_2real_dte, MPIR_2double_dte, 
                    MPIR_2complex_dte, MPIR_2dcomplex_dte, 
                    MPIR_int1_dte, 
                    MPIR_int2_dte, MPIR_int4_dte, MPIR_real4_dte, 
                    MPIR_real8_dte, MPI_REAL, MPI_DOUBLE_PRECISION;
/* 
   The layouts for the types MPI_DOUBLE_INT etc are simply
   struct { 
       double var;
       int    loc;
   }
   
 */
/* Communicators */
typedef struct MPIR_COMMUNICATOR *MPI_Comm;
extern MPI_Comm MPI_COMM_WORLD, MPI_COMM_SELF;

/* Groups */
typedef struct MPIR_GROUP *MPI_Group;
extern MPI_Group MPI_GROUP_EMPTY;

/* Collective operations */
typedef struct MPIR_OP *MPI_Op;
extern MPI_Op MPI_MAX, MPI_MIN, MPI_SUM, MPI_PROD, MPI_LAND, MPI_BAND, 
              MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR, MPI_MINLOC, MPI_MAXLOC;

/* Permanent key values */
/* C Versions (return pointer to value) */
extern int MPI_TAG_UB, MPI_HOST, MPI_IO;
/* Fortran Versions (return value) */
extern int MPIR_TAG_UB, MPIR_HOST, MPIR_IO;

/* Define some null objects */
#define MPI_COMM_NULL      ((MPI_Comm)0)
#define MPI_OP_NULL        ((MPI_Op)0)
#define MPI_GROUP_NULL     ((MPI_Group)0)
#define MPI_DATATYPE_NULL  ((MPI_Datatype)0)
#define MPI_REQUEST_NULL   ((MPI_Request)0)
#define MPI_ERRHANDLER_NULL 0

/* These are only guesses; make sure you change them in mpif.h as well */
#define MPI_MAX_PROCESSOR_NAME 256
#define MPI_MAX_ERROR_STRING   256

/* Pre-defined constants */
#define MPI_UNDEFINED      (-32766)
#define MPI_UNDEFINED_RANK MPI_UNDEFINED
#define MPI_KEYVAL_INVALID 0

/* Upper bound on the overhead in bsend for each message buffer */
#define MPI_BSEND_OVERHEAD 512

/* Topology types */
#define MPI_GRAPH  1
#define MPI_CART   2

#define MPI_BOTTOM      (void *)0

#define MPI_PROC_NULL   (-1)
#define MPI_ANY_SOURCE 	(-2)
#define MPI_ANY_TAG	(-1)


/* Status object.  It is the only user-visible MPI data-structure */
typedef struct { 
    int count;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
} MPI_Status;

/* Must be able to hold any valid address.  64 bit machines may need
   to change this */
typedef long MPI_Aint;


#if (defined(__STDC__) || defined(__cpluscplus))
typedef void (MPI_Handler_function)( MPI_Comm *, int *, ... );
#else
typedef void (MPI_Handler_function)();
#endif

typedef struct MPIR_Errhandler *MPI_Errhandler;
extern MPI_Errhandler MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN, 
       MPIR_ERRORS_WARN;
/* Make the C names for the null functions all upper-case.  Note that 
   this is required for systems that use all uppercas names for Fortran 
   externals.  */
#define MPI_NULL_COPY_FN MPIR_null_copy_fn
#define MPI_NULL_DELETE_FN MPIR_null_delete_fn
#define MPI_DUP_FN MPIR_dup_fn

/* MPI request opjects */
typedef union MPIR_HANDLE *MPI_Request;

/* User combination function */
#if defined(__STDC__) || defined(__cplusplus)
typedef void (MPI_User_function)( void *invec, void *inoutvec, int *len,  
				   MPI_Datatype *datatype); 
#else
typedef void (MPI_User_function)();
#endif

/* MPI Attribute copy and delete functions */
#if defined(__STDC__) || defined(__cplusplus)
typedef int (MPI_Copy_function)( MPI_Comm *oldcomm, int *keyval, 
				 void *extra_state,
			         void *attr_in, void *attr_out, int *flag);
typedef int (MPI_Delete_function)( MPI_Comm *comm, int *keyval, void *attr_val,
			           void *extra_state );
#else
typedef int (MPI_Copy_function)( );
typedef int (MPI_Delete_function)( );
#endif

/* MPI's error classes */
#include "mpi_errno.h"

/* Bindings of the MPI routines */
#include "binding.h"

#if defined(__cplusplus)
};
#endif

#endif




