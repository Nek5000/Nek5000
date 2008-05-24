/*
 * 
 * $Copyright
 * Copyright 1993, 1994, 1995  Intel Corporation
 * INTEL CONFIDENTIAL
 * The technical data and computer software contained herein are subject
 * to the copyright notices; trademarks; and use and disclosure
 * restrictions identified in the file located in /etc/copyright on
 * this system.
 * Copyright$
 * 
 */
 
/*
 * nx.h
 *
 * Copyright (c) 1992 Intel Corporation
 */

#ifndef __NX_H__
#define __NX_H__

#include <standards.h>

_SSD_FUNC_PROTO_BEGIN_


#include <sys/types.h>
#include <sys/estat.h>
#include <pfs/pfs.h>
#include <allocsys.h>

/*
 * Invalid ptype
 */
#define INVALID_PTYPE	-1

/*
 * First number in the Force type range
 */
#define FORCE_TYPE	0x40000000	/* 1073741824 decimal */

/*
 * Message passing info array
 */
extern long msginfo[8];

/*
 * Hardware clock rate
 */
#define HWHZ	10000000

#define HOST		myhost()

#ifndef	_ESIZE_T
#define	_ESIZE_T
/*
 * PFS extended size structure.
 */
struct s_size {
	long	slow;
	long	shigh;
};
typedef struct s_size esize_t;
#endif

#ifndef	_SIZE_SET
#define	_SIZE_SET
/*
 * Symbolic constants for the lsize() and esize() functions.
 */
#define SIZE_SET 0	/* set file size to MAX(cur size,"offset") */
#define SIZE_CUR 1	/* set file size to MAX(cur size,cur offset+"offset" */
#define SIZE_END 2	/* set file size to MAX(cur size,cur size+"offset" */
#endif

#ifndef _IOMODES_SET
#define _IOMODES_SET
/*
 * Symbolic constants for the setiomode() and gopen() functions.
 */
#define M_UNIX   0	/* Unshared file pointer, (default). 
			 */
#define M_LOG    1	/* Shared file pointer, accessed in first come
			 * first served order.
			 */
#define M_SYNC   2	/* Shared file pointer, variable record length, 
			 * accessed in logical node order.
			 */
#define M_RECORD 3 	/* Unshared file pointer, fixed record length,
			 * access in first come first served order, records
			 * are position dependent by logical node number.
			 */
#define M_GLOBAL 4	/* Shared file pointer, all nodes perform the 
			 * same i/o operation. 
			 */
#define M_ASYNC  5	/* Unshared file pointer, no synchronization between
			 * nodes, unrestricted file access
			 */
#endif

#ifndef	NX_ATTR_END
#define NX_ATTR_END	-1
/*
 * Symbolic constants for nx_mkpart_attr and nx_initve_attr
*/
#define	NX_ATTR_SZ	1	/* # of node for partition or application
				 */
#define	NX_ATTR_RECT	2	/* HxW for partition or application
				 */
#define	NX_ATTR_MAP	3	/* List of nodes for partition or application
				 */
#define	NX_ATTR_ANCHOR	4	/* Upper-left corner of rectangle
				 */
#define	NX_ATTR_SCHED	5	/* Partition's scheduling type
				 */
#define	NX_ATTR_RELAXED	6	/* Relaxed allocation flag
				 */
#define	NX_ATTR_EPL	7	/* Effective priority limit
				 */
#define	NX_ATTR_RQ	8	/* Rollin quantum
				 */
#define	NX_ATTR_MOD	9	/* Protection mode
				 */
#define	NX_ATTR_SEL	10	/* Node attribute selection string
				 */
#define	NX_ATTR_PRI	11	/* Priority
				 */
#define	NX_ATTR_PKT	12	/* Packet size
				 */
#define	NX_ATTR_MBF	13	/* Message buffers
				 */
#define	NX_ATTR_MEX	14	/* Memory export
				 */
#define	NX_ATTR_MEA	15	/* Memory each
				 */
#define	NX_ATTR_NOC	16	/* # of correspondents
				 */
#define	NX_ATTR_STH	17	/* Send threshold
				 */
#define	NX_ATTR_SCT	18	/* Send count
				 */
#define	NX_ATTR_GTH	19	/* Give threshold
				 */
#define	NX_ATTR_PLK	20	/* Process lock into memory
				 */
#define	NX_ATTR_ACCT	21	/* MACS Account info
				 */
#define	NX_ATTR_GUEST	22	/* Guest operating system
				 */

/* Parameters for NX_ATTR... attributes
 */
#define NX_SUNMOS	1	/* For NX_ATTR_GUEST, the guest is Sunmos */

#endif


#ifndef _NX_SCHEDTYPE_SET
#define _NX_SCHEDTYPE_SET
/*
 * Symbolic partition scheduling type constants for nx_mkpart(), 
 * nx_mkpart_rect(), and nx_mkpart_map() functions.
 */
#define NX_SPS   2 	/* Use space sharing only for partition  */
#define NX_GANG  1 	/* Use gang scheduling for partition     */
#define NX_STD   0      /* Use standard scheduling for partition */
#endif


extern void	cprobe __protomagic(( long));
extern void	cprobex __protomagic(( long, long, long, long *));
extern void	crecv __protomagic(( long, void *, long));
extern void	crecvx __protomagic(( long, void *, long, long, long, long *));
extern void	csend __protomagic(( long, void *, long, long, long));
extern long	csendrecv __protomagic(( long, void *, long, long, long, long, void *, long));
extern double	dclock __protomagic(( void));
extern void	flick __protomagic(( void));
extern void	flushmsg __protomagic(( long, long, long));
extern void	gcol __protomagic(( void *, long, void *, long, long *));
extern void	gcolx __protomagic(( void *, long *, void *));
extern void	gdhigh __protomagic(( double *, long, double *));
extern void	gdlow __protomagic(( double *, long, double *));
extern void	gdprod __protomagic(( double *, long, double *));
extern void	gdsum __protomagic(( double *, long, double *));
extern void	giand __protomagic(( long *, long, long *));
extern void	gihigh __protomagic(( long *, long, long *));
extern void	gilow __protomagic(( long *, long, long *));
extern void	gior __protomagic(( long *, long, long *));
extern void	giprod __protomagic(( long *, long, long *));
extern void	gisum __protomagic(( long *, long, long *));
extern void	gland __protomagic(( long *, long, long *));
extern void	glor __protomagic(( long *, long, long *));
extern int	gopen __protomagic(( const char *, int, int, mode_t));
extern void	gopf __protomagic((void *, long, void *, long (*function)()));
extern void	gsendx __protomagic(( long, void *, long, long *, long));
extern void	gshigh __protomagic(( float *, long, float *));
extern void	gslow __protomagic(( float *, long, float *));
extern void	gsprod __protomagic(( float *, long, float *));
extern void	gssum __protomagic(( float *, long, float *));
extern void	gsync __protomagic(( void));
extern void	hrecv __protomagic(( long, void *, long, void (*handler)()));
extern void	hrecvx __protomagic(( long, void *, long, long, long, void (* xhandler)(),
			long));
extern void	hsend __protomagic(( long, void *, long, long, long, void (*handler)()));
extern void	hsendrecv __protomagic(( long, void *, long, long, long, long, void *, long,
			   void (*handler)()));
extern void	hsendx __protomagic(( long, void *, long, long, long, void (*xhandler)(),
			long));
extern long	infocount __protomagic(( void));
extern long	infonode __protomagic(( void));
extern long	infoptype __protomagic(( void));
extern long	infotype __protomagic(( void));
extern long	iprobe __protomagic(( long));
extern long	iprobex __protomagic(( long, long, long, long *));
extern long	irecv __protomagic(( long, void *, long));
extern long	irecvx __protomagic(( long, void *, long, long, long, long *));
extern long	isend __protomagic(( long, void *, long, long, long));
extern long	isendrecv __protomagic(( long, void *, long, long, long, long, void *, long));
extern long	masktrap __protomagic(( long));
extern void	msgcancel __protomagic(( long));
extern long	msgdone __protomagic(( long));
extern void	msgignore __protomagic(( long));
extern long	msgmerge __protomagic(( long, long));
extern void	msgwait __protomagic(( long));
extern long	myhost __protomagic(( void));
extern long	mynode __protomagic(( void));
extern long	myptype __protomagic(( void));
extern long	numnodes __protomagic(( void));
extern long	nx_app_nodes __protomagic(( pid_t, nx_nodes_t *, unsigned long *));
extern long	nx_app_rect __protomagic(( long *, long *));
extern long	nx_chpart_epl __protomagic(( char *, long));
extern long	nx_chpart_mod __protomagic(( char *, long));
extern long	nx_chpart_name __protomagic(( char *, char *));
extern long	nx_chpart_owner __protomagic(( char *, long, long));
extern long	nx_chpart_rq __protomagic(( char *, long));
extern long	nx_chpart_sched __protomagic(( char *, long));
extern long	nx_initve __protomagic(( char *, long, char *, int *, char **));
extern long	nx_initve_rect __protomagic(( char *, long, long, long, char *, int *,
				char **));
extern long	nx_initve_attr __protomagic(( char *, int *, char**, ...));
extern long	nx_load __protomagic(( long *, long, long, long *, char *));
extern long	nx_loadve __protomagic(( long *, long, long, long *, char *, char **,
			   char **));
extern long	nx_mkpart __protomagic(( char *, long, long));
extern long	nx_mkpart_attr __protomagic(( char *, ...));
extern long	nx_mkpart_map __protomagic(( char *, long, long *, long));
extern long	nx_mkpart_rect __protomagic(( char *, long, long, long));
extern long	nx_nfork __protomagic(( long *, long, long, long *));
extern void	nx_perror __protomagic(( char *));
extern long	nx_pri __protomagic(( long, long));
extern int	nx_pspart __protomagic(( char *, nx_pspart_t **, unsigned long *));
extern long	nx_rmpart __protomagic(( char *, long, long));
extern long	nx_waitall __protomagic(( void));
extern void	setptype __protomagic(( long));

/* underscore versions */
extern long	_cprobe __protomagic(( long));
extern long	_cprobex __protomagic(( long, long, long, long *));
extern long	_crecv __protomagic(( long, void *, long));
extern long	_crecvx __protomagic(( long, void *, long, long, long, long *));
extern long	_csend __protomagic(( long, void *, long, long, long));
extern long	_csendrecv __protomagic(( long, void *, long, long, long, long, void *, long));
extern double	_dclock __protomagic(( void));
extern long	_flick __protomagic(( void));
extern long	_flushmsg __protomagic(( long, long, long));
extern long	_gcol __protomagic(( void *, long, void *, long, long *));
extern long	_gcolx __protomagic(( void *, long *, void *));
extern long	_gdhigh __protomagic(( double *, long, double *));
extern long	_gdlow __protomagic(( double *, long, double *));
extern long	_gdprod __protomagic(( double *, long, double *));
extern long	_gdsum __protomagic(( double *, long, double *));
extern long	_giand __protomagic(( long *, long, long *));
extern long	_gihigh __protomagic(( long *, long, long *));
extern long	_gilow __protomagic(( long *, long, long *));
extern long	_gior __protomagic(( long *, long, long *));
extern long	_giprod __protomagic(( long *, long, long *));
extern long	_gisum __protomagic(( long *, long, long *));
extern long	_gland __protomagic(( long *, long, long *));
extern long	_glor __protomagic(( long *, long, long *));
extern int	_gopen __protomagic(( const char *, int, int, mode_t));
extern long	_gopf __protomagic((void *, long, void *, long (*function)()));
extern long	_gsendx __protomagic(( long, void *, long, long *, long));
extern long	_gshigh __protomagic(( float *, long, float *));
extern long	_gslow __protomagic(( float *, long, float *));
extern long	_gsprod __protomagic(( float *, long, float *));
extern long	_gssum __protomagic(( float *, long, float *));
extern long	_gsync __protomagic(( void));
extern long	_hrecv __protomagic(( long, void *, long, void (*handler)()));
extern long	_hrecvx __protomagic(( long, void *, long, long, long, void (*xhandler)(),
			long));
extern long	_hsend __protomagic(( long, void *, long, long, long, void (*handler)()));
extern long	_hsendrecv __protomagic(( long, void *, long, long, long, long, void *, long,
			   void (*handler)()));
extern long	_hsendx __protomagic(( long, void *, long, long, long, void (*xhandler)(),
			long));
extern long	_infocount __protomagic(( void));
extern long	_infonode __protomagic(( void));
extern long	_infoptype __protomagic(( void));
extern long	_infotype __protomagic(( void));
extern long	_iprobe __protomagic(( long));
extern long	_iprobex __protomagic(( long, long, long, long *));
extern long	_irecv __protomagic(( long, void *, long));
extern long	_irecvx __protomagic(( long, void *, long, long, long, long *));
extern long	_isend __protomagic(( long, void *, long, long, long));
extern long	_isendrecv __protomagic(( long, void *, long, long, long, long, void *, long));
extern long	_masktrap __protomagic(( long));
extern long	_msgcancel __protomagic(( long));
extern long	_msgdone __protomagic(( long));
extern long	_msgignore __protomagic(( long));
extern long	_msgmerge __protomagic(( long, long));
extern long	_msgwait __protomagic(( long));
extern long	_myhost __protomagic(( void));
extern long	_mynode __protomagic(( void));
extern long	_myptype __protomagic(( void));
extern long	_numnodes __protomagic(( void));
extern long	_setptype __protomagic(( long));

/* iPSC and Touchstone DELTA Compatibility Calls */
extern long	ginv __protomagic(( long));
extern long	gray __protomagic(( long));
extern void	hwclock __protomagic(( esize_t *));
extern long	infopid __protomagic(( void));
extern void	killcube __protomagic(( long, long));
extern void	killproc __protomagic(( long, long));
extern void	led __protomagic(( long));
/*  comment out prototype for load() to avoid conflict with loader.h
extern long	load( char *, long, long); */
extern unsigned long	mclock __protomagic(( void));
extern long	mypart __protomagic(( long *, long *));
extern long	mypid __protomagic(( void));
extern long	nodedim __protomagic(( void));
extern long	restrictvol __protomagic(( int, int, int *));

/* underscore versions */
extern long	_ginv __protomagic(( long));
extern long	_gray __protomagic(( long));
extern void	_hwclock __protomagic(( esize_t *));
extern long	_infopid __protomagic(( void));
extern int	_killcube __protomagic(( long, long));
extern int	_killproc __protomagic(( long, long));
extern long	_led __protomagic(( long));
extern long	_load __protomagic(( char *, long, long));
extern unsigned long	_mclock __protomagic(( void));
extern long	_nodedim __protomagic(( void));
extern long	_restrictvol __protomagic(( int, int, int *));


/*
 * PFS extern definitions.
 */
#ifndef _KERNEL
extern void	cread __protomagic(( int, void *, unsigned int));
extern void	creadv __protomagic(( int, struct iovec *, int));
extern void	cwrite __protomagic(( int, void *, unsigned int));
extern void	cwritev __protomagic(( int, struct iovec *, int));
extern esize_t	eadd __protomagic(( esize_t, esize_t));
extern esize_t	eadd1 __protomagic(( esize_t, long));
extern long	ecmp __protomagic(( esize_t, esize_t));
extern long	ediv __protomagic(( esize_t, long));
extern long	emod __protomagic(( esize_t, long));
extern esize_t	emul __protomagic(( esize_t, long));
extern esize_t	eseek __protomagic(( int, esize_t, int));
extern esize_t	esize __protomagic(( int, esize_t, int));
extern long	estat __protomagic(( char *, struct estat *));
extern esize_t	esub __protomagic(( esize_t, esize_t));
extern esize_t	esub1 __protomagic(( esize_t, long));
extern void	etos __protomagic(( esize_t, char *));
extern long	festat __protomagic(( int, struct estat *));
extern long	fstatpfs __protomagic(( int, struct estatfs *, struct statpfs *,
			  unsigned int));
extern long	getpfsinfo __protomagic(( struct pfsmntinfo **));
extern long	iodone __protomagic(( long));
extern long	iomode __protomagic(( int));
extern void	iowait __protomagic(( long));
extern long	iread __protomagic(( int, void *, unsigned int));
extern long	ireadv __protomagic(( int, struct iovec *, int));
/* ireadoff, ireadvoff newly added	*/
extern long	ireadoff __protomagic(( int, esize_t, char *, unsigned int));
extern long	ireadvoff __protomagic(( int, esize_t, struct iovec *, int));
extern long	iseof __protomagic(( int));
extern long	iwrite __protomagic(( int, void *, unsigned int));
extern long	iwritev __protomagic(( int, struct iovec *, int));
/* iwriteoff, iwritevoff newly added	*/
extern long	iwriteoff __protomagic(( int, esize_t, char *, unsigned int));
extern long	iwritevoff __protomagic(( int, esize_t, struct iovec *, int));
extern long	lestat __protomagic(( char *, struct estat *));
extern long	lsize __protomagic(( int, off_t, int));
/* readoff, readvoff newly added	*/
extern long	readoff __protomagic(( int, esize_t, char *, unsigned int));
extern long	readvoff __protomagic(( int, esize_t, struct iovec *, int));
extern void	setiomode __protomagic(( int, int));
extern long	statpfs __protomagic(( char *, struct estatfs *, struct statpfs *,
			 unsigned int));
extern esize_t	stoe __protomagic(( char *));

/* writeoff, writevoff newly added	*/
extern long	writeoff __protomagic(( int, esize_t, char *, unsigned int));
extern long	writevoff __protomagic(( int, esize_t, struct iovec *, int));

/* underscore versions */
extern int	_cread __protomagic(( int, void *, unsigned int));
extern int	_creadv __protomagic(( int, struct iovec *, int));
extern int	_cwrite __protomagic(( int, void *, unsigned int));
extern int	_cwritev __protomagic(( int, struct iovec *, int));
extern esize_t	_eadd __protomagic(( esize_t, esize_t));
extern esize_t	_eadd1 __protomagic(( esize_t, long));
extern long	_ecmp __protomagic(( esize_t, esize_t));
extern long	_ediv __protomagic(( esize_t, long));
extern long	_emod __protomagic(( esize_t, long));
extern esize_t	_emul __protomagic(( esize_t, long));
extern esize_t	_eseek __protomagic(( int, esize_t, int));
extern esize_t	_esize __protomagic(( int, esize_t, int));
extern long	_estat __protomagic(( char *, struct estat *));
extern esize_t	_esub __protomagic(( esize_t, esize_t));
extern esize_t	_esub1 __protomagic(( esize_t, long));
extern long	_etos __protomagic(( esize_t, char *));
extern long	_festat __protomagic(( int, struct estat *));
extern long	_iodone __protomagic(( long));
extern long	_iomode __protomagic(( int));
extern int	_iowait __protomagic(( long));
extern long	_iread __protomagic(( int, void *, unsigned int));
extern long	_ireadv __protomagic(( int, struct iovec *, int));
extern long	_iseof __protomagic(( int));
extern long	_iwrite __protomagic(( int, void *, unsigned int));
extern long	_iwritev __protomagic(( int, struct iovec *, int));
extern long	_lestat __protomagic(( char *, struct estat *));
extern long	_lsize __protomagic(( int, off_t, int));
extern long	_setiomode __protomagic(( int, int));
extern esize_t	_stoe __protomagic(( char *));

#endif /* _KERNEL */


_SSD_FUNC_PROTO_END_

#endif /* __NX_H__ */

