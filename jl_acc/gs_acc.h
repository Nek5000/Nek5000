#ifndef GS_ACC_H
#define GS_ACC_H

#ifndef GS_ACC_C
extern void gs_flatmap_setup_acc(const sint *handle, int n, struct gs_data **fgs_info);

extern void fgs_fields_acc(const sint *handle, double *u, const sint *stride, const sint *n,
			   const sint *dom, const sint *op, const sint *transpose,
			   struct gs_data **fgs_info);
#endif

#endif
