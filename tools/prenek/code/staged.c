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
