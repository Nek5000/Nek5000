/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
sgl_iadd(int *vals, int level)
{
  int edge, type, dest, source, len, mask, ceil;
  long msg_id;
  int tmp, *work;


  /* all msgs will be of the same length */
  work = &tmp;
  len = INT_LEN;

  if (level > i_log2_num_nodes)
    {error_msg_fatal("sgl_add() :: level too big?");}

  if (level<=0)
    {return;}

  /* implement the mesh fan in/out exchange algorithm */
  if (my_id<floor_num_nodes)
    {
      mask = 0;
      for (edge = 0; edge < level; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = edge_node[edge];
	      type = 10001 + my_id + (num_nodes*edge);
	      if (my_id > dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  vals[0] += work[0];
		}
	    }
	  mask <<= 1;
	  mask += 1;
	}
    }

  if (my_id<floor_num_nodes)
    {
      mask >>= 1;
      for (edge = 0; edge < level; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = edge_node[level-edge-1];
	      type = 10001 + my_id + (num_nodes*edge);
	      if (my_id < dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  vals[0] = work[0];
		}
	    }
	  mask >>= 1;
	}
    }
}



/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
******************************************************************************/
void
sgl_radd(double *vals, int level)
{
  int edge, type, dest, source, len, mask, ceil;
  long msg_id;
  double tmp, *work;


  /* all msgs will be of the same length */
  work = &tmp;
  len = FLOAT_LEN;

  if (level > i_log2_num_nodes)
    {error_msg_fatal("sgl_add() :: level too big?");}

  if (level<=0)
    {return;}

  /* implement the mesh fan in/out exchange algorithm */
  if (my_id<floor_num_nodes)
    {
      mask = 0;
      for (edge = 0; edge < level; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = edge_node[edge];
	      type = 1000001 + my_id + (num_nodes*edge);
	      if (my_id > dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  vals[0] += work[0];
		}
	    }
	  mask <<= 1;
	  mask += 1;
	}
    }

  if (my_id<floor_num_nodes)
    {
      mask >>= 1;
      for (edge = 0; edge < level; edge++)
	{
	  if (!(my_id & mask))
	    {
	      source = dest = edge_node[level-edge-1];
	      type = 1000001 + my_id + (num_nodes*edge);
	      if (my_id < dest)
		{
		  msg_id = isend(type,vals,len,dest,0);
		  msgwait(msg_id);
		}
	      else
		{
		  type =  type - my_id + source;
		  msg_id = irecv(type,work,len);
		  msgwait(msg_id);
		  vals[0] = work[0];
		}
	    }
	  mask >>= 1;
	}
    }
}  





/******************************************************************************
Function: giop()

Input : 
Output: 
Return: 
Description: 
 
ii+1 entries in seg :: 0 .. ii

******************************************************************************/
void 
ssgl_radd(register double *vals, register double *work, register int level, 
	  register int *segs)
{
  register int edge, type, dest, source, mask;
  register int stage_n;

#ifdef DEBUG
  if (level > i_log2_num_nodes)
    {error_msg_fatal("sgl_add() :: level > log_2(P)!!!");}
#endif

  /* all msgs are *NOT* the same length */
  /* implement the mesh fan in/out exchange algorithm */
  for (mask=0, edge=0; edge<level; edge++, mask++)
    {
      stage_n = (segs[level] - segs[edge]);
      if (stage_n && !(my_id & mask))
	{
	  dest = edge_node[edge];
	  type = 100001 + my_id + (num_nodes*edge);
	  if (my_id>dest)
	    {csend(type, vals+segs[edge],stage_n<<3,dest,0);}
	  else
	    {
	      type =  type - my_id + dest;
	      crecv(type,work,stage_n<<3);
	      rvec_add(vals+segs[edge], work, stage_n); 
/*            daxpy(vals+segs[edge], work, stage_n); */
	    }
	}
      mask <<= 1;
    }
  mask>>=1;
  for (edge=0; edge<level; edge++)
    {
      stage_n = (segs[level] - segs[level-1-edge]);
      if (stage_n && !(my_id & mask))
	{
	  dest = edge_node[level-edge-1];
	  type = 10000001 + my_id + (num_nodes*edge);
	  if (my_id<dest)
	    {csend(type,vals+segs[level-1-edge],stage_n<<3,dest,0);}
	  else
	    {
	      type =  type - my_id + dest;
	      crecv(type,vals+segs[level-1-edge],stage_n<<3);
	    }
	}
      mask >>= 1;
    }
}  


/***********************************comm.c*************************************
Function: grop()

Input : 
Output: 
Return: 
Description: fan-in/out version

note good only for num_nodes=2^k!!!

***********************************comm.c*************************************/
void
grop_hc_vvl(REAL *vals, REAL *work, int *segs, int *oprs, int dim)
{
  register int mask, edge, n;
  int type, dest;
  vfp fp;
#if defined MPISRC
  MPI_Status  status;
#elif defined NXSRC
  int len;
#endif


#ifdef SAFE
  /* ok ... should have some data, work, and operator(s) */
  if (!vals||!work||!oprs||!segs)
    {error_msg_fatal("grop_hc() :: vals=%d, work=%d, oprs=%d",vals,work,oprs);}

  /* non-uniform should have at least two entries */
  if ((oprs[0] == NON_UNIFORM)&&(n<2))
    {error_msg_fatal("grop_hc() :: non_uniform and n=0,1?");}    

  /* check to make sure comm package has been initialized */
  if (!p_init)
    {comm_init();}
#endif

  /* if there's nothing to do return */
  if ((num_nodes<2)||(dim<=0))
    {return;}

  /* the error msg says it all!!! */
  if (modfl_num_nodes)
    {error_msg_fatal("grop_hc() :: num_nodes not a power of 2!?!");}

  /* can't do more dimensions then exist */
  dim = MIN(dim,i_log2_num_nodes);

  /* advance to list of n operations for custom */
  if ((type=oprs[0])==NON_UNIFORM)
    {oprs++;}

  if ((fp = (vfp) rvec_fct_addr(type)) == NULL)
    {
      error_msg_warning("grop_hc() :: hope you passed in a rbfp!\n");
      fp = (vfp) oprs;
    }

#if  defined NXSRC
  for (mask=1,edge=0; edge<dim; edge++,mask<<=1)
    {
      n = segs[dim]-segs[edge];
      len = n*REAL_LEN;
      dest = my_id^mask;
      if (my_id > dest)
	{csend(76207+my_id,(char *)vals,len,dest,0); break;}
      else
	{crecv(76207+dest,(char *)work,len); (*fp)(vals, work, n, oprs);}
    }
      
  if (edge==dim)
    {mask>>=1;}
  else
    {while (++edge<dim) {mask<<=1;}}
  
  for (edge=0; edge<dim; edge++,mask>>=1)
    {
      if (my_id%mask)
	{continue;}
      len = (segs[dim]-segs[dim-1-edge])*REAL_LEN;
      
      dest = my_id^mask;
      if (my_id < dest)
	{csend(163841+my_id,(char *)vals,len,dest,0);}
      else
	{crecv(163841+dest,(char *)vals,len);}
    }

#elif defined MPISRC 
  for (mask=1,edge=0; edge<dim; edge++,mask<<=1)
    {
      n = segs[dim]-segs[edge];
      dest = my_id^mask;
      if (my_id > dest)
	{MPI_Send(vals,n,REAL_TYPE,dest,76207+my_id,MPI_COMM_WORLD);}
      else
	{
	  MPI_Recv(work,n,REAL_TYPE,MPI_ANY_SOURCE,76207+dest,
		   MPI_COMM_WORLD, &status);
	  (*fp)(vals, work, n, oprs);
	}
    }

  if (edge==dim)
    {mask>>=1;}
  else
    {while (++edge<dim) {mask<<=1;}}

  for (edge=0; edge<dim; edge++,mask>>=1)
    {
      if (my_id%mask)
	{continue;}
      
      n = (segs[dim]-segs[dim-1-edge]);
      
      dest = my_id^mask;
      if (my_id < dest)
	{MPI_Send(vals,n,REAL_TYPE,dest,163841+my_id,MPI_COMM_WORLD);}
      else
	{
	  MPI_Recv(vals,n,REAL_TYPE,MPI_ANY_SOURCE,163841+dest,
		   MPI_COMM_WORLD, &status);
	}
    }
#else
  return;
#endif
}  


/* look in xxt_model comm.c for source */

/***********************************comm.c*************************************
Function: grop_hc_ ()

Input : 
Output: 
Return: 
Description: fan-in/out version (Fortran)
***********************************comm.c*************************************/
#if defined UPCASE
void
GROP_HC_VVL  (REAL *vals, REAL *work, int *len_vec, int *oprs, int *dim)
#else
void
grop_hc_vvl_ (REAL *vals, REAL *work, int *len_vec, int *oprs, int *dim)
#endif
{
  grop_hc_vvl(vals, work, len_vec, oprs, *dim);
}
