C
C  This file is part of MUMPS 5.4.1, released
C  on Tue Aug  3 09:49:43 UTC 2021
C
C
C  Copyright 1991-2021 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
      MODULE MUMPS_STATIC_MAPPING
      USE MUMPS_LR_COMMON
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: MUMPS_DISTRIBUTE, MUMPS_RETURN_CANDIDATES,
     &     MUMPS_INIT_ARCH_PARAMETERS,MUMPS_END_ARCH_CV
      integer,pointer,dimension(:,:),SAVE::cv_cand
      integer,pointer,dimension(:),SAVE::cv_par2_nodes
      integer,SAVE::cv_slavef,cv_nb_niv2,cv_lp,cv_mp
      integer, parameter:: tsplit_beg=4 
      integer, parameter:: tsplit_mid=5 
      integer, parameter:: tsplit_last=6 
      integer,parameter::cv_invalid=-9999
      DOUBLE PRECISION,parameter::cv_d_invalid=-9999.D0
      integer,parameter::cv_equilib_flops=1
      integer,parameter::cv_equilib_mem=2
      integer,parameter::cv_error_memalloc = -13
      integer,parameter::cv_error_memdeloc = -96
      integer,dimension(:),allocatable,save :: mem_distribtmp
      integer, dimension(:),allocatable, save :: table_of_process
      integer,dimension(:),allocatable,save :: mem_distribmpi
      integer, save ::ke69,nb_arch_nodes
      logical,dimension(:),allocatable,save :: allowed_nodes
      integer,dimension(:),allocatable,save :: score
      type nodelist
         integer::nodenumber
         type(nodelist),pointer::next
      end type nodelist
      type alloc_arraytype
         integer, pointer, dimension(:)::t2_nodenumbers
         integer, pointer, dimension(:,:)::t2_cand
         DOUBLE PRECISION, pointer, dimension(:)::t2_candcostw(:),
     &                                            t2_candcostm(:)
         integer:: nmb_t2s
      end type alloc_arraytype
      type splitting_data
         integer:: new_ison,new_ifather,old_keep2
         DOUBLE PRECISION:: ncostw_oldinode,ncostm_oldinode,
     &                      tcostw_oldinode,tcostm_oldinode
      end type splitting_data
      type procs4node_t
         integer, dimension(:), pointer :: ind_proc
      end type procs4node_t
      DOUBLE PRECISION, pointer, dimension(:) ::
     &      cv_proc_workload,
     &      cv_proc_maxwork,
     &      cv_proc_memused,
     &      cv_proc_maxmem
      type(splitting_data)::cv_last_splitting
      integer::cv_n,cv_nsteps,cv_maxlayer,
     &         cv_nbsa,cv_maxnsteps,cv_maxdepth,
     &         cv_maxnodenmb,cv_total_amalg,cv_total_split,
     &         cv_bitsize_of_int,cv_size_ind_proc
     &         ,cv_mixed_strat_bound,cv_dist_L0_mixed_strat_bound
     &         ,cv_layerl0_end,cv_layerl0_start
      integer :: layerL0_endforarrangeL0
      DOUBLE PRECISION :: mincostw
      DOUBLE PRECISION:: cv_costw_upper,cv_costm_upper,
     &      cv_costw_layer0,cv_costm_layer0,cv_relax,
     &      cv_costw_total,cv_costm_total,cv_l0wthresh,cv_splitthresh
      logical::cv_constr_work,cv_constr_mem
      integer,pointer,dimension(:):: cv_nodetype,cv_nodelayer,
     &          cv_layerl0_array,cv_proc_sorted,cv_depth
      integer,dimension(:),pointer::
     &    cv_ne,cv_nfsiz,cv_frere,cv_fils,cv_keep,cv_info,
     &    cv_procnode,cv_ssarbr,cv_icntl
      integer(8),dimension(:),pointer::cv_keep8
      type(alloc_arraytype),pointer,dimension(:)::cv_layer_p2node
      DOUBLE PRECISION,dimension(:),pointer:: cv_ncostw,
     & cv_tcostw,cv_ncostm,cv_tcostm,cv_layerworkload,cv_layermemused
     & ,cv_layerl0_sorted_costw
      type(procs4node_t),dimension(:),pointer :: cv_prop_map
      integer, dimension(:), pointer :: cv_SIZEOFBLOCKS
      logical  :: cv_BLKON
      contains
      subroutine MUMPS_DISTRIBUTE(n,slavef,icntl,info,
     &                        ne,nfsiz,frere,fils,keep,KEEP8,
     &                        procnode,ssarbr,nbsa,peak,istat
     &                        , SIZEOFBLOCKS, LSIZEOFBLOCKS
     &     )
      implicit none
      integer,intent(in)::n,slavef
      integer, intent(inout),TARGET:: ne(n),nfsiz(n),
     &         procnode(n),ssarbr(n),frere(n),fils(n),keep(500),
     &         icntl(60),info(80)
      integer, intent(in) :: LSIZEOFBLOCKS
      integer, intent(in) :: SIZEOFBLOCKS(LSIZEOFBLOCKS)
      INTEGER(8) KEEP8(150)
      integer,intent(out)::nbsa,istat
      integer ierr,nmb_thislayer,layernmb,mapalgo,allocok,i
      integer,pointer,dimension(:)::thislayer
      integer,parameter::memonly=1,floponly=2,hybrid=3
      DOUBLE PRECISION::
     &         maxwork,minwork,maxmem,minmem,workbalance,membalance
      DOUBLE PRECISION:: cost_root_node
      DOUBLE PRECISION,dimension(:),allocatable:: work_per_proc
      integer,dimension(:),allocatable::id_son
      logical :: cont
      character (len=48):: err_rep,subname
      DOUBLE PRECISION peak
      logical :: BLKON
      BLKON    = (SIZEOFBLOCKS(1).GT.0)
      cv_BLKON = BLKON
      istat=-1
      subname='DISTRIBUTE'
      cv_lp=icntl(1)
      cv_mp=icntl(3)
      IF (icntl(4).LT.2) cv_mp=0
      nullify(thislayer)
      err_rep='INITPART1'
      call MUMPS_INITPART1(n,slavef,
     &                   frere,fils,nfsiz,ne,keep,KEEP8,icntl,info,
     &                   procnode,ssarbr,peak,ierr
     &                      , SIZEOFBLOCKS, LSIZEOFBLOCKS
     &     )
      if (ierr.ne.0) goto 99999
      err_rep='PROCINIT'
      call MUMPS_PROCINIT(istat=ierr)
      if (ierr.ne.0) goto 99999
      err_rep='CALCCOST'
      call MUMPS_CALCCOSTS(ierr)
      if (ierr.ne.0) goto 99999
      err_rep='ROOTLIST'
      call MUMPS_ROOTLIST(ierr)
      if (ierr.ne.0) goto 99999
      err_rep='LAYERL0'
      call MUMPS_LAYERL0(ierr)
      if (ierr.ne.0) goto 99999
      if (ierr.ne.0) goto 99999
      err_rep='INITPART2'
      call MUMPS_INITPART2(ierr)
      if (ierr.ne.0) goto 99999
      err_rep='WORKMEM_'
      call MUMPS_WORKMEM_IMBALANCE(
     &     cv_proc_workload,cv_proc_memused,
     &     maxwork,minwork,maxmem,minmem)
      if(maxwork.gt.0.0D0) then
         workbalance=minwork/maxwork
      else
         workbalance=0.0D0
      endif
      if(maxmem.gt.0.0D0) then
         membalance=minmem/maxmem
      else
         membalance=0.0D0
      endif
      err_rep='mem_alloc'
      allocate(thislayer(cv_maxnodenmb),STAT=allocok)
      if (allocok.gt.0) then
         cv_info(1) = cv_error_memalloc
         cv_info(2) = 2*cv_maxnsteps+cv_maxnodenmb
         if(cv_lp.gt.0)
     &   write(cv_lp,*)'memory allocation error in ',subname
         ierr = cv_error_memalloc
         goto 99999
      end if
      cont=.TRUE. 
      layernmb=0
      mapalgo=floponly  
      err_rep='SELECT_TYPE3'
      call MUMPS_SELECT_TYPE3(ierr)
      if (ierr.ne.0) goto 99999
      IF (cv_keep(38) .ne. 0 .and. cv_keep(60) .eq. 0 ) THEN
        call MUMPS_GET_FLOPS_COST(cv_nfsiz(keep(38)),
     &               cv_nfsiz(keep(38)), cv_nfsiz(keep(38)),
     &               cv_keep(50), 3, cost_root_node)
        cost_root_node = cost_root_node / dble(cv_slavef)
        do i=1, cv_slavef
          cv_proc_memused(i)=cv_proc_memused(i)+
     &        dble(cv_nfsiz(keep(38)))*dble(cv_nfsiz(keep(38)))/
     &        dble(cv_slavef)
          cv_proc_workload(i)=cv_proc_workload(i)+dble(cost_root_node)
        enddo
      ENDIF
      do while((cont).OR.(layernmb.le.cv_maxlayer))
         err_rep='FIND_THIS'
         call MUMPS_FIND_THISLAYER(layernmb,thislayer,nmb_thislayer,
     &                                ierr)
         if (ierr.ne.0) goto 99999
         err_rep='DO_SPLITTING'
         if(cv_keep(82) .gt. 0) then
            if(layernmb.gt.0) call MUMPS_SPLIT_DURING_MAPPING
     &           (layernmb,thislayer,nmb_thislayer,ierr)
         endif
         if (ierr.ne.0) goto 99999
         err_rep='ASSIGN_TYPES'
         call MUMPS_ASSIGN_TYPES(layernmb,thislayer,nmb_thislayer,
     &                              ierr)
         if (ierr.ne.0) goto 99999
         if(layernmb.gt.0) then
            if ((cv_keep(24).eq.1).OR.(cv_keep(24).eq.2).OR.
     &          (cv_keep(24).eq.4).OR.(cv_keep(24).eq.6)) then
               err_rep='COSTS_LAYER_T2'
               call MUMPS_COSTS_LAYER_T2(layernmb,nmb_thislayer,ierr)
            elseif((cv_keep(24).eq.8).OR.(cv_keep(24).eq.10)
     &             .OR.(cv_keep(24).eq.12).OR.(cv_keep(24).eq.14)
     &             .OR.(cv_keep(24).eq.16).OR.(cv_keep(24).eq.18)) then
               err_rep='COSTS_LAYER_T2PM'
               call MUMPS_COSTS_LAYER_T2PM(layernmb,nmb_thislayer,ierr)
            else
               err_rep='wrong strategy for COSTS_LAYER_T2'
               ierr = -9999
            endif
            if (ierr.ne.0) goto 99999
            err_rep='WORKMEM_'
            call MUMPS_WORKMEM_IMBALANCE(
     &                          cv_proc_workload,cv_proc_memused,
     &                                   maxwork,minwork,maxmem,minmem)
            if(maxwork.gt.0.0D0) then
               workbalance=minwork/maxwork
            else
               workbalance=0.0D0
            endif
            if(maxmem.gt.0.0D0) then
               membalance=minmem/maxmem
            else
               membalance=0.0D0
            endif
            if(mapalgo.eq.memonly) then
               err_rep='MAP_LAYER'
               call MUMPS_MAP_LAYER(layernmb,thislayer,
     &              nmb_thislayer,cv_equilib_mem,ierr)
               if (ierr.ne.0) goto 99999
            elseif(mapalgo.eq.floponly) then
               err_rep='MAP_LAYER'
               call MUMPS_MAP_LAYER(layernmb,thislayer,
     &              nmb_thislayer,cv_equilib_flops,ierr)
               if (ierr.ne.0) goto 99999
            elseif(mapalgo.eq.hybrid) then
               if (workbalance <= membalance) then
                  err_rep='MAP_LAYER'
                  call MUMPS_MAP_LAYER(layernmb,thislayer,
     &                 nmb_thislayer,cv_equilib_flops,ierr)
                  if (ierr.ne.0) goto 99999
               else
                  err_rep='MAP_LAYER'
                  call MUMPS_MAP_LAYER(layernmb,thislayer,
     &                 nmb_thislayer,cv_equilib_mem,ierr)
                  if (ierr.ne.0) goto 99999
               endif
            else
               if(cv_lp.gt.0)
     &         write(cv_lp,*)'Unknown mapalgo in ',subname
               return
            endif
      endif
         layernmb=layernmb+1
         err_rep='HIGHER_LAYER'
         call MUMPS_HIGHER_LAYER(layernmb,thislayer,
     &                 nmb_thislayer,cont,ierr)
         if (ierr.ne.0) goto 99999
      end do
      IF ( (cv_keep(79).EQ.0).OR.(cv_keep(79).EQ.3).OR.
     &     (cv_keep(79).EQ.5).OR.(cv_keep(79).EQ.7)
     &   ) THEN
       if(cv_slavef.gt.4) then
          err_rep='POSTPROCESS'
          call MUMPS_POSTPROCESS_MEM()
       endif
      ENDIF 
      err_rep='SETUP_CAND'
      call MUMPS_SETUP_CAND(ierr)
      if (ierr.ne.0) goto 99999
      err_rep='ENCODE_PROC'
      call MUMPS_ENCODE_PROCNODE(ierr)
      if (ierr.ne.0) goto 99999
      err_rep='STORE_GLOB'
      call MUMPS_STORE_GLOBALS(ne,nfsiz,frere,fils,keep,KEEP8,
     &                         info,procnode,ssarbr,nbsa)
      err_rep='mem_dealloc'
      deallocate(thislayer,STAT=ierr)
      if (ierr.ne.0) then
         if(cv_lp.gt.0)
     &   write(cv_lp,*)'Memory deallocation error in ',subname
         ierr = cv_error_memdeloc
         goto 99999
      endif
      err_rep='TERMGLOB'
      call MUMPS_TERMGLOB(ierr)
      if (ierr.ne.0) goto 99999
      istat=0
      return
99999 continue
      if(cv_lp.gt.0) then
         write(cv_lp,*)'Error in ',subname,', layernmb=',layernmb
         write(cv_lp,*)'procedure reporting the error: ',err_rep
      endif
      if(ierr.eq.cv_error_memalloc) then
         info(1) = cv_info(1)
         info(2) = cv_info(2)
      endif
      istat=ierr
      return
      CONTAINS
      subroutine MUMPS_ACCEPT_L0(
     &                       map_strat,workload,memused,accepted,
     &                              istat)
      implicit none
      integer,intent(in)::map_strat
      DOUBLE PRECISION,dimension(:),intent(in)::workload, memused
      logical,intent(out)::accepted
      integer,intent(out)::istat
      DOUBLE PRECISION maxi,mini,mean,stddev, dpkeep102
      integer i,nmb
      intrinsic maxval,minval,count,sum
      character (len=48):: subname
      logical alternative_criterion
      DOUBLE PRECISION::
     &          MINFLOPS , MINMEM,
     &          CL_RATE, DV_RATE
      istat=-1
      if ( cv_keep(72) .EQ. 1) then
       MINFLOPS = 2.0D0
       MINMEM=50.0D0
       CL_RATE =0.8D0
       DV_RATE=0.2D0
      else
       IF (cv_keep(66).NE.0) THEN
         MINFLOPS = 5.0D8
         MINMEM=5.0D7
         CL_RATE =0.8D0
         DV_RATE=0.2D0
       ELSE
         MINFLOPS = 5.0D7
         MINMEM=5.0D6
         CL_RATE =0.8D0
         DV_RATE=0.2D0
       ENDIF
      endif
      dpkeep102 = dble(cv_keep(102))
      IF (cv_keep(66).NE.0) THEN
        IF (cv_slavef.LT.3)THEN
         dpkeep102 = dble(150)
        ELSEIF (cv_slavef.LT.5)THEN
         dpkeep102 = dble(200)
        ELSEIF (cv_slavef.LT.8)THEN
         dpkeep102 = dble(250)
        ELSEIF (cv_slavef.LT.32)THEN
         dpkeep102 = dble(275)
        ELSEIF (cv_slavef.LT.512)THEN
         dpkeep102 = dble(300)
        ELSEIF (cv_slavef.GE.512)THEN
         dpkeep102 = dble(400)
        ENDIF
      ENDIF
      subname='ACCEPT_L0'
      accepted=.FALSE.
      alternative_criterion=.FALSE. 
      if(map_strat.eq.cv_equilib_flops) then
         maxi=maxval(workload)
         mini=minval(workload)
         if (maxi.lt.MINFLOPS) then
            accepted=.TRUE.
         elseif(maxi.le.(dpkeep102/dble(100))*mini)then
            accepted=.TRUE.
         endif
         if ((.NOT.accepted).AND.(alternative_criterion)) then
            mean=sum(workload)/max(dble(cv_slavef),dble(1))
            stddev=dble(0)
            do i=1,cv_slavef
               stddev=stddev+
     &               (abs(workload(i)-mean)*abs(workload(i)-mean))
            enddo
            stddev=sqrt(stddev/max(dble(cv_slavef),dble(1)))
            nmb=count(mask=abs(workload-mean)<stddev)
            if((dble(nmb)/max(dble(cv_slavef),dble(1)).gt.CL_RATE)
     &       .AND.(stddev.lt.DV_RATE*mean)) accepted=.TRUE.
         endif
      elseif(map_strat.eq.cv_equilib_mem) then
         maxi=maxval(memused)
         mini=minval(memused)
         if (maxi.lt.MINMEM) then
            accepted=.TRUE.
         else if(cv_slavef.lt.48) then
            if (maxi.le.dble(2)*mini) accepted=.TRUE.
         else if(cv_slavef.lt.128) then
            if (maxi.le.dble(4)*mini) accepted=.TRUE.
         else if(cv_slavef.lt.256) then
            if (maxi.le.dble(6)*mini) accepted=.TRUE.
         else if(cv_slavef.lt.512) then
            if (maxi.le.dble(8)*mini) accepted=.TRUE.
         else if(cv_slavef.gt.512) then
            if (maxi.le.dble(10)*mini) accepted=.TRUE.
         end if
      endif
      istat=0
      return
      end subroutine MUMPS_ACCEPT_L0
      subroutine MUMPS_ARRANGEL0(map_strat,layerL0end,workload,memused,
     &                           procnode,istat,respect_prop)
      implicit none
      integer, intent(in)::map_strat, layerL0end
      DOUBLE PRECISION,dimension(:),intent(out)::workload, memused
      integer, intent(out)::procnode(:),istat
      logical, intent(in), OPTIONAL:: respect_prop
      integer i,j,ierr, nodenumber,proc
      DOUBLE PRECISION work,mem
      character (len=48):: err_rep,subname
      istat=-1
      subname='ARRANGEL0'
      if ((.NOT.associated(cv_tcostw)).OR.(.NOT.associated(cv_tcostm)))
     &   then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error:tcost must be allocated in ',subname
         return
      end if
      if((map_strat.ne.cv_equilib_flops).and.
     &   (map_strat.ne.cv_equilib_mem)) return
      do i=1,cv_n
         procnode(i)=cv_invalid
      end do
      do i=1,cv_slavef
         workload(i)=cv_proc_workload(i)
         memused(i)=cv_proc_memused(i)
      end do
      do i=cv_layerl0_start,layerL0end   
         nodenumber=cv_layerl0_array(i)
         work=cv_tcostw(nodenumber)
         mem=cv_tcostm(nodenumber)
         err_rep='FIND_BEST_PROC'
         if(present(respect_prop)) then
            call MUMPS_FIND_BEST_PROC(nodenumber,map_strat,work,mem,
     &           workload,memused,proc,ierr,respect_prop)
         else
            call MUMPS_FIND_BEST_PROC(nodenumber,map_strat,work,mem,
     &           workload,memused,proc,ierr)
         endif
         if(ierr.eq.0) then
            procnode(nodenumber)=proc
         else
         if(cv_lp.gt.0)
     &   write(cv_lp,*)'Error reported by ',err_rep,' in ',subname
            do j=1,cv_slavef
               workload(j)=cv_proc_workload(j)
               memused(j)=cv_proc_memused(j)
            end do
            do j=1,cv_n
               procnode(j)=cv_invalid
            end do
            return
         end if
      end do
      istat=0
      return
      end subroutine MUMPS_ARRANGEL0
      subroutine MUMPS_ASSIGN_TYPES( layernmb,thislayer,nmb_thislayer,
     &                                 istat )
      implicit none
      integer,intent(in)::layernmb,thislayer(:), 
     &                    nmb_thislayer
      integer,intent(out)::istat
      integer i,in,npiv,nfront,inode,inoderoot,par_nodes_in_layer,
     &        dummy,allocok
      character (len=48):: subname
      istat=-1
      subname='ASSIGN_TYPES'
      if((layernmb.lt.0).or.(layernmb.gt.cv_maxlayer)) return
      if(cv_slavef.eq.1) then
         if(layernmb.eq.0) then
            do inode=1,cv_n
               cv_nodetype(inode)=0
            end do
         end if
      else if(layernmb.eq.0) then
         do i=1,nmb_thislayer
            inode=thislayer(i)
            inoderoot=inode
            if(cv_nodetype(inode).ne.cv_invalid) cycle
            cv_nodetype(inode)=0
 30         continue
            in = inode
            do while (in .ne. 0)
               inode = in
               do while (in .gt. 0)
                  in = cv_fils(in)
               end do
               if (in.lt.0) in=-in
            end do
 10         continue
            if ( inode .ne. inoderoot ) then
               cv_nodetype(inode)=-1
               in = cv_frere(inode)
               inode = abs(in)
               if (in .lt. 0) then
                  go to 10
               else
                  go to 30
               end if
            end if
         end do
      else
         do i=1,nmb_thislayer
            inode=thislayer(i)
            in = inode
            npiv = 0
            do while (in.gt.0)
               if (cv_BLKON) then
                 npiv = npiv + cv_SIZEOFBLOCKS(in)
               else
                 npiv = npiv + 1
               endif
               in = cv_fils(in)
            end do
            nfront = cv_nfsiz(inode)
            if(cv_nodetype(inode).ne.cv_invalid) cycle 
            if( ( MUMPS_ISTYPE2BYSIZE(nfront,npiv)) .AND.
     &           (in.ne.0)) then
               cv_nodetype(inode)=2
            else
               cv_nodetype(inode)=1
            end if
         end do
      end if
      if(layernmb.gt.0) then
         par_nodes_in_layer=0
         do i=1,nmb_thislayer
            inode=thislayer(i)
            if (MUMPS_IS_NODE_OF_TYPE2(inode))
     &         par_nodes_in_layer=par_nodes_in_layer+1
         enddo
         if(par_nodes_in_layer.gt.0) then
            allocate(
     &cv_layer_p2node(layernmb)%t2_nodenumbers(par_nodes_in_layer),
     &cv_layer_p2node(layernmb)%t2_cand(par_nodes_in_layer,cv_slavef+1),
     &cv_layer_p2node(layernmb)%t2_candcostw(par_nodes_in_layer),
     &cv_layer_p2node(layernmb)%t2_candcostm(par_nodes_in_layer),
     &               STAT=allocok)
            if (allocok.gt.0) then
               cv_info(1) = cv_error_memalloc
               cv_info(2) = (3+cv_slavef+1)*par_nodes_in_layer
               istat = cv_error_memalloc
               if(cv_lp.gt.0)
     &         write(cv_lp,*)'memory allocation error in ',subname
               return
            end if
            cv_layer_p2node(layernmb)%nmb_t2s=par_nodes_in_layer
            dummy=1
            do i=1,nmb_thislayer
               inode=thislayer(i)
               if (MUMPS_IS_NODE_OF_TYPE2(inode)) then
                  cv_layer_p2node(layernmb)%t2_nodenumbers(dummy)=inode
                  cv_layer_p2node(layernmb)%t2_cand(dummy,:)=0
                  cv_layer_p2node(layernmb)%t2_candcostw(dummy)
     &                                                    =cv_d_invalid
                  cv_layer_p2node(layernmb)%t2_candcostm(dummy)
     &                                                    =cv_d_invalid
                  dummy=dummy+1
               endif
            enddo
         else
            nullify(cv_layer_p2node(layernmb)%t2_nodenumbers,
     &              cv_layer_p2node(layernmb)%t2_cand,
     &              cv_layer_p2node(layernmb)%t2_candcostw,
     &              cv_layer_p2node(layernmb)%t2_candcostm)
         end if
      endif
      istat=0
      return
      end subroutine MUMPS_ASSIGN_TYPES
      function MUMPS_BIT_GET(procs4node,procnumber)
      implicit none
      integer,intent(in)::procs4node(cv_size_ind_proc)
      integer,intent(in)::procnumber
      logical :: MUMPS_BIT_GET
      integer pos1,pos2
      pos1 = (procnumber-1)/cv_bitsize_of_int +1
      pos2 = mod(procnumber-1,cv_bitsize_of_int)
      MUMPS_BIT_GET=btest(procs4node(pos1),pos2)
      return
      end function MUMPS_BIT_GET
      function MUMPS_BIT_GET4PROC(inode,procnumber)
!DEC$ NOOPTIMIZE
      implicit none
      integer, intent(in)::inode,procnumber
      logical :: MUMPS_BIT_GET4PROC
      integer pos1,pos2
      MUMPS_BIT_GET4PROC=.FALSE.
      if((procnumber.lt.1).or.(procnumber.gt.cv_slavef)) return
      if(.not.associated(cv_prop_map(inode)%ind_proc)) return
      pos1 = (procnumber-1)/cv_bitsize_of_int +1
      pos2 = mod(procnumber-1,cv_bitsize_of_int)
      MUMPS_BIT_GET4PROC=btest
     &               (cv_prop_map(inode)%ind_proc(pos1),pos2)
      return
      end function MUMPS_BIT_GET4PROC
      subroutine MUMPS_BIT_SET(procs4node,procnumber,istat)
      implicit none
      integer, intent(inout)::procs4node(cv_size_ind_proc)
      integer,intent(in)::procnumber
      integer, intent(out)::istat
      integer pos1,pos2
      istat = -1
      if((procnumber.lt.1).or.(procnumber.gt.cv_slavef)) return
      if(cv_bitsize_of_int.le.0) return
      pos1 = (procnumber-1)/cv_bitsize_of_int +1
      pos2 = mod(procnumber-1,cv_bitsize_of_int)
      procs4node(pos1)=ibset(procs4node(pos1),pos2)
      istat = 0
      return
      end subroutine MUMPS_BIT_SET
      subroutine MUMPS_CALCCOSTS(istat)
      implicit none
      integer,intent(out)::istat
      integer i
      DOUBLE PRECISION :: maxcostw_root
      istat = -1
      if ((.NOT.associated(cv_tcostw)).OR.(.NOT.associated(cv_tcostm)))
     &   then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)
     &        'Error: tcost must be allocated in MUMPS_CALCCOSTS'
         return
      end if
      maxcostw_root = 0D0
      do i=1,cv_n
         if (cv_frere(i).eq.cv_n+1) then
            cv_tcostw(i)=0.0D0
            cv_ncostw(i)=0.0D0
            cv_tcostm(i)=0.0D0
            cv_ncostm(i)=0.0D0
         elseif (cv_frere(i).eq.0) then
            cv_depth(i)=1
            call MUMPS_TREECOSTS(i)
            maxcostw_root = max(maxcostw_root,cv_tcostw(i))
         end if
      end do
      istat = 0
      mincostw = 1.0D0+maxcostw_root/(dble(cv_maxnsteps)*
     &              dble(10*cv_slavef) )
      return
      end subroutine MUMPS_CALCCOSTS
      subroutine MUMPS_CALCNODECOSTS(npiv,nfront,costw,costm)
      implicit none
      integer,intent(in)::npiv,nfront
      DOUBLE PRECISION,intent(out)::costw,costm
      character (len=48):: subname
      subname='CALCNODECOSTS'
      if((npiv.le.1).and.(nfront.le.1)) then
         costw = dble(0)
         costm = dble(1)
      else
         if((cv_keep(494).ne.0).and.(cv_keep(471).ge.0).and.
     &     (npiv.ge.cv_keep(490)).and.(nfront.ge.cv_keep(491))) then
           WRITE(*,*) " *** Temp internal error in MUMPS_CALCNODECOSTS:"
           CALL MUMPS_ABORT()
           call MUMPS_CALCNODECOSTS_BLR(npiv, nfront, costw, costm,
     &              cv_keep(471), cv_keep(472), cv_keep(475), 
     &              cv_keep(488), cv_keep(50))
         else
           if(cv_keep(50).eq.0) then
             costw= 2.0D0*dble(nfront)*dble(npiv)*dble(nfront-npiv-1)
     &      + dble(npiv)*dble(npiv+1)*dble(2*npiv+1)/dble(3)
     &           + dble(2*nfront-npiv-1) * dble(npiv) / dble(2)
             costm= dble(npiv)*(dble(2*nfront)-dble(npiv))
           else
             costw= dble(npiv) *
     &            (dble(nfront)*dble(nfront)+dble(2*nfront) -
     &             dble(nfront+1) * dble(npiv+1) +
     &             dble(npiv+1) * dble(2*npiv+1) / dble(6))
             costm= dble(npiv) * dble(nfront)
           end if
         end if
      end if
      if((costw.lt.0.0D0).or.(costm.lt.0.0D0)) then
      endif
      return
      end subroutine MUMPS_CALCNODECOSTS
      SUBROUTINE MUMPS_CALCNODECOSTS_BLR(NPIV, NFRONT, COSTW, COSTM,
     &                   K471, K472, K475, K488, SYM)
        INTEGER, INTENT(IN) :: NPIV, NFRONT, SYM, K471, K472, K475, K488
        DOUBLE PRECISION, INTENT(OUT) :: COSTW, COSTM
        INTEGER :: IBCKSZ
        DOUBLE PRECISION :: B,R,M,N
        M = DBLE(NPIV)
        N = DBLE(NFRONT)
        CALL COMPUTE_BLR_VCS(K472, IBCKSZ, K488, NPIV) 
        B = DBLE(IBCKSZ)
        B = MIN(B,M)
        IF (K471.EQ.0) THEN
          R = 1.0D0
        ELSEIF (K471.EQ.1) THEN
          R = SQRT(DBLE(N))
       ELSE
          WRITE(*,*) 'Internal error in MUMPS_CALCNODECOSTS_BLR', K471
          CALL MUMPS_ABORT()
        ENDIF
        R = MIN(R,B/2)
        IF (SYM.EQ.0) THEN
          COSTW = M/B * B*(B+1.0D0)*(2.0D0*B+1.0D0)/3.0D0
          IF (K475.EQ.0) THEN
            COSTW = COSTW + 2.0D0*M/(B*B)*(N-(M+B)/2.0D0) * B*B*B
          ELSEIF (K475.EQ.1) THEN
            COSTW = COSTW + M/(B*B)*(N-(M+B)/2.0D0) * B*B*(R+B)
          ELSEIF (K475.EQ.2) THEN
            COSTW = COSTW + M/(B*B)*(2.0D0*N-3.0D0*M-2.0D0*B) * B*B*R
     &                    + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 * B*B*B
          ELSEIF (K475.EQ.3) THEN
            COSTW = COSTW + 2.0D0*M/(B*B)*(N-(M+B)/2.0D0) * B*B*R
          ENDIF
          COSTW = COSTW + 2.0D0*M/(B*B)*(N-(M+B)/2.0D0) * 2.0D0*B*B*R
          COSTW = COSTW + (4.0D0*B*R*R + 2.0D0*B*B*R) * (
     &                  (N-M)*(N-M)*M/(B*B*B) 
     &                + (N-M)/B*(M/B-1.0D0)*M/B 
     &                + (M/B-1.0D0)*M/B*(2.0D0*M/B-1.0D0)/6.0D0 
     &       )     
          COSTM = M*(2.0D0*N-M)/(B*B) * 2.0D0*B*R
        ELSE
          COSTW = M/B * B*(B+1.0D0)*(2.0D0*B+1.0D0)/6.0D0
          IF (K475.EQ.0.OR.K475.EQ.1) THEN
            COSTW = COSTW + M/(B*B)*(N-(M+B)/2.0D0) * B*B*B
          ELSEIF (K475.EQ.2) THEN
            COSTW = COSTW + (N-M)*M/(B*B) * B*B*R
     &                    + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 * B*B*B
          ELSEIF (K475.EQ.3) THEN
            COSTW = COSTW + M/(B*B)*(N-(M+B)/2.0D0) * B*B*R
          ENDIF
          COSTW = COSTW + M/(B*B)*(N-(M+B)/2.0D0) * 2.0D0*B*B*R
          COSTW = COSTW + (4.0D0*B*R*R + 2.0D0*B*B*R) * (
     &                      (N-M)*(N-M)*M/(B*B*B)/2.0D0 
     &                    + (N-M)/B*(M/B-1.0D0)*M/B/2.0D0 
     &                    + (M/B-1.0D0)*M/B*(M/B+1.0D0)/6.0D0 
     &       )     
          COSTM = M*N/(B*B) * 2.0D0*B*R
        ENDIF
      END SUBROUTINE MUMPS_CALCNODECOSTS_BLR
      SUBROUTINE MUMPS_COSTS_BLR_T2_MASTER(NPIV, NFRONT, 
     &                   COSTW, COSTM, K471, K472, K475, K488, SYM)
        INTEGER, INTENT(IN) :: NPIV, NFRONT, SYM, K471, K472, K475, K488
        DOUBLE PRECISION, INTENT(OUT) :: COSTW, COSTM
        INTEGER :: IBCKSZ
        DOUBLE PRECISION :: B,R,M,N
        M = DBLE(NPIV)
        N = DBLE(NFRONT)
        CALL COMPUTE_BLR_VCS(K472, IBCKSZ, K488, NPIV) 
        B = DBLE(IBCKSZ)
        B = MIN(B,M)
        IF (K471.EQ.0) THEN
          R = 1.0D0
        ELSEIF (K471.EQ.1) THEN
          R = SQRT(DBLE(N))
       ELSE
          WRITE(*,*) 'Internal error in ',
     &       'MUMPS_COSTS_BLR_T2_MASTER', K471
          CALL MUMPS_ABORT()
        ENDIF
        R = MIN(R,B/2)
        IF (SYM.EQ.0) THEN
          COSTW = M/B * B*(B+1.0D0)*(2.0D0*B+1.0D0)/3.0D0
          IF (K475.EQ.0) THEN
            COSTW = COSTW + M/(B*B)*(N-(M+B)/2.0D0) * B*B*B 
     &                + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 * B*B*B 
          ELSEIF (K475.EQ.1) THEN
            COSTW = COSTW + M/(B*B)*(N-(M+B)/2.0D0) * B*B*B 
     &                + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 * B*B*R 
          ELSEIF (K475.EQ.2) THEN
            COSTW = COSTW + M/(B*B)*(N-M) * B*B*R 
     &                + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 * B*B*(B+R)
          ELSEIF (K475.EQ.3) THEN
            COSTW = COSTW + M/(B*B)*(N-(M+B)/2.0D0) * B*B*R 
     &                + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 * B*B*R 
          ENDIF
          COSTW = COSTW + M/(B*B)*(N-(M+B)/2.0D0) * 2.0D0*B*B*R 
     &              + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 * 2.0D0*B*B*R 
          COSTW = COSTW + (4.0D0*B*R*R + 2.0D0*B*B*R) * (
     &                   (N-M)/B*(M/B-1.0D0)*M/B/2.0D0 
     &                 + (M/B-1.0D0)*M/B*(2.0D0*M/B-1.0D0)/6.0D0 
     &       )     
          COSTM = M*N/(B*B) * 2.0D0*B*R
        ELSE
          COSTW = M/B * B*(B+1.0D0)*(2.0D0*B+1.0D0)/6.0D0
          IF (K475.LE.2) THEN
            COSTW = COSTW + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 * B*B*B
          ELSEIF (K475.EQ.3) THEN
            COSTW = COSTW + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 * B*B*R
          ENDIF
          COSTW = COSTW + (M/B-1.0D0)*M/B*(M/B-1.0D0)/6.0D0 
     &                       * 2.0D0*B*B*R 
          COSTW = COSTW + (4.0D0*B*R*R + 2.0D0*B*B*R) * (
     &                   (M/B-1.0D0)*M/B*(M/B+1.0D0)/6.0D0 
     &       )     
          COSTM = M*M/(B*B) * 2.0D0*B*R
        ENDIF
      END SUBROUTINE MUMPS_COSTS_BLR_T2_MASTER
      SUBROUTINE MUMPS_COSTS_BLR_T2_SLAVE(NPIV, NFRONT, 
     &                NROW, COSTW, COSTM, K471, K472, K475, K488, SYM)
        INTEGER, INTENT(IN) :: NPIV, NFRONT, SYM, K471, K472, 
     &                         K475, K488
        DOUBLE PRECISION, INTENT(IN) :: NROW
        DOUBLE PRECISION, INTENT(OUT) :: COSTW, COSTM
        INTEGER :: IBCKSZ
        DOUBLE PRECISION :: B,R,M,N,P
        M = NROW
        N = DBLE(NFRONT)
        P = DBLE(NPIV)
        CALL COMPUTE_BLR_VCS(K472, IBCKSZ, K488, NPIV) 
        B = DBLE(IBCKSZ)
        B = MIN(B,M)
        IF (K471.EQ.0) THEN
          R = 1.0D0
        ELSEIF (K471.EQ.1) THEN
          R = SQRT(DBLE(N))
       ELSE
          WRITE(*,*) 'Internal error in ',
     &       'MUMPS_COSTS_BLR_T2_SLAVE', K471
          CALL MUMPS_ABORT()
        ENDIF
        R = MIN(R,B/2)
        COSTW = 0.0D0
        IF (K475.EQ.0) THEN
          COSTW = COSTW + (M*P)/(B*B) * B*B*B 
        ELSE 
          COSTW = COSTW + (M*P)/(B*B) * B*B*R 
        ENDIF
        COSTW = COSTW + (M*P)/(B*B) * 2.0D0*B*B*R
        IF (SYM.EQ.0) THEN
          COSTW = COSTW + (4.0D0*B*R*R + 2.0D0*B*B*R) * (
     &                  M/B*(P/B-1.0D0)*P/B/2.0D0 
     &                + (N-M)*M*P/(B*B*B) 
     &       )     
        ELSE
          COSTW = COSTW + (4.0D0*B*R*R + 2.0D0*B*B*R) * (
     &                  M/B*(P/B-1.0D0)*P/B/2.0D0 
     &                + (N-M)*M*P/(B*B*B)/2.0D0 
     &       )     
        ENDIF
        COSTM = M*P/(B*B) * 2.0D0*B*R
      END SUBROUTINE MUMPS_COSTS_BLR_T2_SLAVE
      subroutine MUMPS_COSTS_LAYER_T2(layernmb,nmb_thislayer,istat)
      implicit none
      integer,intent(in)::layernmb,nmb_thislayer
      integer,intent(out)::istat
      integer in,inode,j,kmax,npiv,nfront,ncb,ncol,
     &        min_needed,max_needed,more_than_needed,total_nmb_cand,
     &        nmb_type2_thislayer,fraction,
     &        total_cand_layer,cand_strat, keep48_loc
      DOUBLE PRECISION flop1,work_type2_thislayer,
     &        relative_weight,workmaster,nrow
      logical force_cand
      intrinsic count,max
      character (len=48):: subname
      integer MUMPS_REG_GETKMAX, MUMPS_BLOC2_GET_NSLAVESMIN,
     &        MUMPS_BLOC2_GET_NSLAVESMAX
      external MUMPS_REG_GETKMAX, MUMPS_BLOC2_GET_NSLAVESMIN,
     &        MUMPS_BLOC2_GET_NSLAVESMAX
      istat=-1
      subname='COSTS_LAYER_T2'
      if (cv_keep(24).lt.1) then
          if(cv_lp.gt.0)
     &    write(cv_lp,*)'Error in ',subname,'. Wrong keep24'
         return
      endif
      force_cand=(mod(cv_keep(24),2).eq.0)
      cand_strat=cv_keep(24)/2
      nmb_type2_thislayer=cv_layer_p2node(layernmb)%nmb_t2s
      if (nmb_type2_thislayer.gt.0) then
         work_type2_thislayer=0.0D0
         do j=1,nmb_type2_thislayer
            inode=cv_layer_p2node(layernmb)%t2_nodenumbers(j)
            work_type2_thislayer=work_type2_thislayer+cv_ncostw(inode)
         end do
         if(cv_relax.le.0.0D0) then
            if(cv_lp.gt.0)
     &           write(cv_lp,*)'Error in ',subname,'. Wrong cv_relax'
            return
         endif
         total_cand_layer=int(cv_relax*dble(cv_slavef))
         do j=1,nmb_type2_thislayer
            inode=cv_layer_p2node(layernmb)%t2_nodenumbers(j)
            nfront=cv_nfsiz(inode)
            npiv=0
            in=inode
            do while(in.gt.0)
               if (cv_BLKON) then
                 npiv = npiv + cv_SIZEOFBLOCKS(in)
               else
                 npiv=npiv+1
               endif
               in=cv_fils(in)
            end do
            ncb=nfront-npiv
            kmax = MUMPS_REG_GETKMAX(cv_keep8(21),ncb)
            if (force_cand) then
              if (cv_keep(50) ==  0) then
                keep48_loc=0
              else
                keep48_loc=3
              endif
              if (cv_keep(48).EQ.5) keep48_loc = 5
               min_needed = MUMPS_BLOC2_GET_NSLAVESMIN(
     &             cv_slavef, keep48_loc,cv_keep8(21),
     &             cv_keep(50),nfront,ncb,
     &             cv_keep(375), cv_keep(119))
               max_needed = MUMPS_BLOC2_GET_NSLAVESMAX(
     &             cv_slavef, keep48_loc,cv_keep8(21),
     &             cv_keep(50),nfront,ncb,
     &             cv_keep(375), cv_keep(119))
               if(cand_strat.eq.1) then
                  more_than_needed = 0
               elseif (cand_strat.eq.2) then
                  if(work_type2_thislayer.gt.0.0D0) then
                  relative_weight=cv_ncostw(inode)/work_type2_thislayer
                  else
                     relative_weight = 0.0D0
                  endif
                  fraction=nint(relative_weight *
     &                 dble(total_cand_layer))
                  more_than_needed=min(max(0,cv_slavef-1-min_needed),
     &                                 max(0,fraction-min_needed)    )
               elseif (cand_strat.eq.3) then
                  more_than_needed=cv_slavef-1-min_needed
               else
                  if(cv_lp.gt.0)
     &            write(cv_lp,*)'Unknown cand. strategy in ',subname
                  return
               endif
               total_nmb_cand=min(min_needed+more_than_needed,
     &              cv_slavef-1)
               total_nmb_cand=min(total_nmb_cand,max_needed)
            else
               total_nmb_cand=0
            endif
            cv_layer_p2node(layernmb)%t2_cand(j,cv_slavef+1)
     &                                 = total_nmb_cand
            if(cv_keep(50).eq.0) then
               flop1=dble(2*npiv)*dble(nfront)-
     &              dble(npiv+nfront)*dble(npiv+1)
               flop1= dble(npiv)*flop1 +
     &          dble(2 * npiv-npiv-1)*dble(npiv)/dble(2)+
     &          dble(npiv)*dble(npiv+1)*dble(2*npiv+1)/dble(3)
            else
               flop1=dble(npiv)*
     &         ( dble(npiv)*dble(npiv)+dble(npiv)-
     &         dble(npiv*npiv+npiv+1) )+
     &         (dble(npiv)*dble(npiv+1)*dble(2*npiv+1))/dble(6)
            endif
            cv_ncostw(inode)=flop1
            if(total_nmb_cand.gt.0) then
               nrow = dble(max(min(dble(ncb)/dble(total_nmb_cand),
     &                     dble(kmax)),
     &                 dble(ncb)/dble(cv_slavef-1)))
            elseif(cv_slavef.gt.1) then
               nrow = dble(max(dble(kmax),
     &                         dble(ncb)/dble(cv_slavef-1)))
            else
               nrow = dble(ncb)
            endif
            if(cv_keep(50).eq.0) then
               flop1 = dble(npiv)*dble(nrow)+
     &                 dble(nrow)*dble(npiv)*dble(2*nfront-npiv-1)
            else
               ncol= nfront   
               flop1 = dble(npiv)*dble(nrow)*
     &          (dble(2*ncol)-dble(nrow)-dble(npiv)+dble(1))
               workmaster = dble(npiv)*dble(npiv)*dble(npiv)/dble(3)
               if (workmaster.gt.flop1) flop1=workmaster
            endif
            cv_layer_p2node(layernmb)%t2_candcostw(j)=flop1
            if(cv_keep(50).eq.0) then
               cv_ncostm(inode)=dble(npiv)*dble(nfront)
            else
               cv_ncostm(inode)=dble(npiv)*dble(npiv)
            endif
            if(cv_keep(50).eq.0) then
               cv_layer_p2node(layernmb)%t2_candcostm(j)
     &                                  =dble(npiv)*dble(nrow)
            else
               cv_layer_p2node(layernmb)%t2_candcostm(j)
     &                                  =dble(npiv)*dble(nrow)
            endif
         end do
      endif
      istat=0
      return
      end subroutine MUMPS_COSTS_LAYER_T2
      subroutine MUMPS_COSTS_LAYER_T2PM(layernmb,nmb_thislayer,istat)
!DEC$ OPTIMIZE:1
      implicit none
      integer,intent(in)::layernmb,nmb_thislayer
      integer,intent(out)::istat
      integer in,inode,j,jj,kmax,npiv,nfront,ncb,ncol,
     &        total_nmb_cand,nmb_type2_thislayer,
     &        total_cand_layer,npropmap,min_needed,
     &        keep48_loc
      DOUBLE PRECISION flop1,work_type2_thislayer,
     &        relative_weight,workmaster,nrow
      DOUBLE PRECISION save_ncostw, save_ncostm
      LOGICAL SPLITNODE, BLRNODE
      intrinsic count,max
      character (len=48):: subname
      integer MUMPS_REG_GETKMAX, MUMPS_BLOC2_GET_NSLAVESMIN
      external MUMPS_REG_GETKMAX, MUMPS_BLOC2_GET_NSLAVESMIN
      istat=-1
      SPLITNODE=.FALSE.
      BLRNODE=.FALSE.
      save_ncostw = 1.0D0
      save_ncostm = 1.0D0
      subname='COSTS_LAYER_T2PM'
      if((cv_keep(24).ne.8).AND.(cv_keep(24).ne.10)
     &    .AND.(cv_keep(24).ne.12).AND.(cv_keep(24).ne.14)
     &    .AND.(cv_keep(24).ne.16).AND.(cv_keep(24).ne.18)) then
          if(cv_lp.gt.0)
     &    write(cv_lp,*)'Error in ',subname,'. Wrong keep24'
         return
      endif
      nmb_type2_thislayer=cv_layer_p2node(layernmb)%nmb_t2s
      if (nmb_type2_thislayer.gt.0) then
         total_cand_layer=0
         work_type2_thislayer=0.0D0
         do j=1,nmb_type2_thislayer
            inode=cv_layer_p2node(layernmb)%t2_nodenumbers(j)
            work_type2_thislayer=work_type2_thislayer+cv_ncostw(inode)
            npropmap=0
            do jj=1,cv_slavef
               if( MUMPS_BIT_GET4PROC(inode,jj))
     &              npropmap=npropmap+1
            end do
            total_cand_layer=total_cand_layer+npropmap
         end do
         do j=1,nmb_type2_thislayer
            inode=cv_layer_p2node(layernmb)%t2_nodenumbers(j)
            nfront=cv_nfsiz(inode)
            SPLITNODE = (abs(cv_nodetype(inode)).GT.3) 
            IF (SPLITNODE) THEN
              save_ncostw = cv_ncostw(inode)
              save_ncostm = cv_ncostm(inode)
            ENDIF
            npiv=0
            in=inode
            do while(in.gt.0)
               if (cv_BLKON) then
                 npiv = npiv + cv_SIZEOFBLOCKS(in)
               else
                 npiv=npiv+1
               endif
               in=cv_fils(in)
            end do
            ncb=nfront-npiv
            kmax = MUMPS_REG_GETKMAX(cv_keep8(21),ncb)
            if(kmax.lt.1) then
               kmax = max(kmax,1)
            endif
            if (cv_keep(50) ==  0) then
              keep48_loc=0
            else
              keep48_loc=3
            endif
            if (cv_keep(48).EQ.5) keep48_loc = 5
            min_needed= MUMPS_BLOC2_GET_NSLAVESMIN
     &          (cv_slavef, keep48_loc,cv_keep8(21),
     &           cv_keep(50),nfront,ncb,
     &           cv_keep(375), cv_keep(119))
            if(min_needed.lt.1) then
               if(cv_lp.gt.0)
     &         write(cv_lp,*)'Error in ',subname,'.NEG min_needed'
               return
            endif
            if ((cv_keep(24).eq.8).OR.(cv_keep(24).eq.14).OR.
     &          (cv_keep(24).eq.18)) then
               npropmap=0
               do jj=1,cv_slavef
                  if( MUMPS_BIT_GET4PROC(inode,jj))
     &                 npropmap=npropmap+1
               end do
               total_nmb_cand=max(npropmap-1,min_needed)
            elseif(cv_keep(24).eq.10) then
               if(work_type2_thislayer.gt.0.0D0) then
                  relative_weight=cv_ncostw(inode)/work_type2_thislayer
               else
                  relative_weight = 0.0D0
               endif
               total_nmb_cand=nint(relative_weight *
     &              dble(total_cand_layer))
               total_nmb_cand=max(total_nmb_cand-1,min_needed)
            elseif((cv_keep(24).eq.12).OR.(cv_keep(24).eq.16)) then
               if(layernmb.lt.cv_dist_L0_mixed_strat_bound) then
                  if(cv_mp.gt.0)then
                     write(cv_mp,*)'Strat', cv_keep(24),
     &                             ': use 8 on layer',layernmb
                  endif
                  npropmap=0
                  do jj=1,cv_slavef
                     if( MUMPS_BIT_GET4PROC(inode,jj))
     &                    npropmap=npropmap+1
                  end do
                  total_nmb_cand=max(npropmap-1,min_needed)
               else
                  if(cv_mp.gt.0)then
                     write(cv_mp,*)'Strat', cv_keep(24),
     &                             ': use 10 on layer',layernmb
                  endif
                  if(work_type2_thislayer.gt.0.0D0) then
                  relative_weight=cv_ncostw(inode)/work_type2_thislayer
                  else
                     relative_weight = 0.0D0
                  endif
                  total_nmb_cand=nint(relative_weight *
     &                 dble(total_cand_layer))
                  total_nmb_cand=max(total_nmb_cand-1,min_needed)
               endif
            else
               if(cv_lp.gt.0)
     &              write(cv_lp,*)'Unknown cand. strategy in ',subname
               return
            endif
            total_nmb_cand=max(total_nmb_cand,1)
            total_nmb_cand=min(total_nmb_cand,cv_slavef-1)
            total_nmb_cand=min(total_nmb_cand,ncb)
            cv_layer_p2node(layernmb)%t2_cand(j,cv_slavef+1)
     &                                 = total_nmb_cand
            BLRNODE = ((cv_keep(494).ne.0).and.(cv_keep(471).ge.0) .and.
     &     (npiv.ge.cv_keep(490)).and.(nfront.ge.cv_keep(491))) 
            IF (BLRNODE) THEN
              call MUMPS_COSTS_BLR_T2_MASTER(npiv, nfront, 
     &              cv_ncostw(inode), cv_ncostm(inode), 
     &              cv_keep(471), cv_keep(472), cv_keep(475), 
     &              cv_keep(488), cv_keep(50))
            ELSE
              if(cv_keep(50).eq.0) then
                flop1=dble(2*npiv)*dble(nfront)-
     &              dble(npiv+nfront)*dble(npiv+1)
                flop1= dble(npiv)*flop1 +
     &          dble(2 * npiv-npiv-1)*dble(npiv)/dble(2)+
     &          dble(npiv)*dble(npiv+1)*dble(2*npiv+1)/dble(3)
              else
                flop1=dble(npiv)*
     &         ( dble(npiv)*dble(npiv)+dble(npiv)-
     &         dble(npiv*npiv+npiv+1) )+
     &         (dble(npiv)*dble(npiv+1)*dble(2*npiv+1))/dble(6)
              endif
              cv_ncostw(inode)=flop1
            ENDIF
            IF (SPLITNODE) THEN
             cv_layer_p2node(layernmb)%t2_candcostw(j)= 
     &       max(save_ncostw - cv_ncostw(inode), 1.0D0)
            ELSE
              if(total_nmb_cand.gt.0) then
                nrow = dble(max(min(dble(ncb)/dble(total_nmb_cand),
     &                     dble(kmax)),
     &                 dble(ncb)/dble(cv_slavef-1)))
              elseif(cv_slavef.gt.1) then
                nrow = dble(max(dble(kmax),
     &                         dble(ncb)/dble(cv_slavef-1)))
              else
                nrow = dble(ncb)
              endif
              IF (BLRNODE) THEN
                call MUMPS_COSTS_BLR_T2_SLAVE(npiv, nfront, 
     &              nrow,           
     &              cv_layer_p2node(layernmb)%t2_candcostw(j),
     &              cv_layer_p2node(layernmb)%t2_candcostm(j),
     &              cv_keep(471), cv_keep(472), cv_keep(475), 
     &              cv_keep(488), cv_keep(50))
              ELSE
                if(cv_keep(50).eq.0) then
                  flop1 = dble(npiv)*dble(nrow)+
     &                 dble(nrow)*dble(npiv)*dble(2*nfront-npiv-1)
                else
                  ncol= nfront   
                  flop1 = dble(npiv)*dble(nrow)*
     &          (dble(2*ncol)-dble(nrow)-dble(npiv)+dble(1))
                  workmaster = dble(npiv)*dble(npiv)*dble(npiv)/dble(3)
                  if (workmaster.gt.flop1) flop1=workmaster
                endif
                cv_layer_p2node(layernmb)%t2_candcostw(j)=flop1
              ENDIF
            ENDIF
            IF (.NOT.BLRNODE) THEN
              if(cv_keep(50).eq.0) then
                cv_ncostm(inode)=dble(npiv)*dble(nfront)
              else
                cv_ncostm(inode)=dble(npiv)*dble(npiv)
              endif
            ENDIF
            IF (SPLITNODE) THEN
              cv_layer_p2node(layernmb)%t2_candcostm(j) = 
     &          max(save_ncostm - cv_ncostm(inode), 1.0D0)
            ELSEIF (.NOT.BLRNODE) THEN
              if(cv_keep(50).eq.0) then
                cv_layer_p2node(layernmb)%t2_candcostm(j)
     &                                  =dble(npiv)*dble(nrow)
              else
                cv_layer_p2node(layernmb)%t2_candcostm(j)
     &                                  =dble(npiv)*dble(nrow)
              endif
            ENDIF
          end do
        endif
        istat=0
        return
      end subroutine MUMPS_COSTS_LAYER_T2PM
      subroutine MUMPS_SPLIT_DURING_MAPPING(
     &     layernmb,thislayer,nmb_thislayer,
     &     istat )
      implicit none
      integer,intent(in)::layernmb,nmb_thislayer
      integer,intent(in)::thislayer(:) 
      integer,intent(out)::istat
      integer i,j,k1,k2,k3,ierr,inode,nfront,npiv,
     &     npropmap, inode_tmp, allocok
      logical doit
      integer, allocatable, dimension(:) :: npivsplit
      integer :: lnpivsplit
      integer :: bsize
      integer :: k1_temp, npiv_beg, npiv_end
      character (len=48):: err_rep,subname
      istat=-1
      subname='SPLIT_DURING_MAPPING'
      if((layernmb.lt.0).or.(layernmb.gt.cv_maxlayer)) return
      if (cv_slavef.eq.1) then
         return
      endif
      if (cv_icntl(59) .ne. 0) then
         istat = 0
         return
      endif
      lnpivsplit = cv_keep(108)
      allocate(npivsplit(lnpivsplit),stat=allocok)
      if (allocok .NE. 0) then
        cv_info(1) = cv_error_memalloc
        cv_info(2) = lnpivsplit
        istat = cv_error_memalloc
        if(cv_lp.gt.0)
     &       write(cv_lp,*)'memory allocation error in ',subname
        return
      endif
      do i=1,nmb_thislayer
        ierr=0
        inode=thislayer(i)
        nfront = cv_nfsiz(inode)
        inode_tmp=inode
        npiv=0
        do while (inode_tmp.gt.0)
          if (cv_BLKON) then
            npiv = npiv + cv_SIZEOFBLOCKS(inode_tmp)
          else
            npiv=npiv+1
          endif
           inode_tmp=cv_fils(inode_tmp)
        end do
        if (inode_tmp .eq. 0) cycle
        npropmap=0
        do j=1,cv_slavef
           if( MUMPS_BIT_GET4PROC(inode,j)) then
              npropmap=npropmap+1
           endif
        end do
        IF ((keep(376) .EQ.1)
     &     ) THEN
           err_rep='GET_SPLIT_4_PERF'
           CALL MUMPS_GET_SPLIT_4_PERF(inode, nfront, npiv,
     &     dble(npropmap),
     &     k1, lnpivsplit, npivsplit, n, cv_frere(1), 
     &     cv_keep(1),
     &     cv_fils(1), cv_BLKON, cv_SIZEOFBLOCKS(1),
     &     istat)
           k3=k1
           doit = .true.
           GOTO 200
        ENDIF
        IF ((cv_keep(79) .EQ.0).OR.(cv_keep(79).GE.5)) THEN
          err_rep='GET_SPLIT_INKPART'
         call MUMPS_GET_SPLIT_INKPART(inode,
     &        doit,npiv,nfront,npropmap,k1,k3,
     &        ierr)
        ELSE
         err_rep='GET_MEMSPLIT_INKPART'
         call MUMPS_GET_MEMSPLIT_INKPART(inode,
     &        doit,npiv,nfront,npropmap,k2,ierr)
         k1=k2
         k3=k2
        ENDIF
        if (ierr.eq.0) then
           if (lnpivsplit < k1) then
            write(*,*) 'error in', subname, lnpivsplit, k1, cv_keep(108)
            call MUMPS_ABORT()
          endif
          bsize   = max(npiv/k1,1)
          if (cv_BLKON) then
            inode_tmp = inode
            npiv_beg = 0  
            npiv_end = 0
            k1_temp  = 0
            do while (inode_tmp.gt.0)
             npiv_end = npiv_end + cv_SIZEOFBLOCKS(inode_tmp)
             if (npiv_end-npiv_beg.ge.bsize) then
               k1_temp = k1_temp+1
               npivsplit(k1_temp) = npiv_end-npiv_beg
               npiv_beg = npiv_end
               if ( ( (npiv-npiv_beg).gt.0) .and.
     &              (npiv-npiv_beg.LT.2*bsize) 
     &            ) then
                k1_temp = k1_temp+1
                npivsplit(k1_temp) = npiv - npiv_beg
                exit
               endif
             endif
             inode_tmp=cv_fils(inode_tmp)
            enddo
            if (k1_temp.eq.0) then
              k1_temp = 1
              npivsplit(1) = npiv
            else
             if (npiv_end.gt.npiv_beg) then
              k1_temp = k1_temp+1
              npivsplit(k1_temp) = npiv_end-npiv_beg
             endif
            endif
            k1 = k1_temp
          else
            do j = 1, k1-1
              npivsplit(j)= bsize
            enddo
            npivsplit(k1) = npiv-bsize*(k1-1)
          endif
        endif
 200    CONTINUE
        if(ierr.ne.0) then
            if(cv_lp.gt.0)
     &           write(cv_lp,*)'Error reported by ',
     &           err_rep,' in ',subname
            istat =ierr
            goto 100
        endif
         if ( ( k1.le.1).or.(k3.le.1).or.(.NOT.doit) ) cycle
         err_rep='SPLITNODE_INKPART'
         call MUMPS_SPLITNODE_INTREE( inode, nfront, npiv, k1,
     &        lnpivsplit, npivsplit, cv_keep(1), n, cv_fils(1),
     &        cv_frere(1),
     &        cv_nfsiz(1), cv_ne(1), cv_info(5), 
     &        cv_nsteps, cv_nodetype(1), ierr
     &                      , SIZEOFBLOCKS, LSIZEOFBLOCKS
     &                      , BLKON
     &      )
         if(ierr.ne.0) then
            if(cv_lp.gt.0)
     &      write(cv_lp,*)'Error reported by ',err_rep,
     &           ' in ',subname
            istat = ierr
            goto 100
         endif
         err_rep='SPLITNODE_UPDATE'
         call MUMPS_SPLITNODE_UPDATE( inode, nfront, npiv, k1,
     &        lnpivsplit, npivsplit,
     &        ierr)
         if(ierr.ne.0) then
            if(cv_lp.gt.0)
     &      write(cv_lp,*)'Error reported by ',err_rep,
     &           ' in ',subname
            istat = ierr
            goto 100
         endif
      end do
      istat=0
 100  continue
      deallocate(npivsplit)
      return
      end subroutine MUMPS_SPLIT_DURING_MAPPING
      subroutine MUMPS_GET_SPLIT_INKPART(inode,
     &     doit,npiv,nfront,npropmap,k1,k3,istat)
      implicit none
      integer,intent(in)::inode
      logical,intent(out)::doit
      integer,intent(in) :: npiv, nfront, npropmap
      integer,intent(out) :: istat
      integer,intent(out) ::k1,k3
      integer npiv2,nfront2,npiv_son2
      integer ncb,kmax,keep48_loc,nslaves_max,
     &     nslaves_estim,strat,kk
      DOUBLE PRECISION wk_master,wk_master2,wk_slave2
      integer MUMPS_REG_GETKMAX,
     &     MUMPS_BLOC2_GET_NSLAVESMAX,
     &     MUMPS_BLOC2_GET_NSLAVESMIN
      external MUMPS_REG_GETKMAX
      external MUMPS_BLOC2_GET_NSLAVESMAX
      external MUMPS_BLOC2_GET_NSLAVESMIN
      doit=.FALSE.
      k1=1
      k3 =1
      istat=-1
      doit=.TRUE.
      if (cv_nodetype(inode) .gt. 0) then
         doit=.FALSE.
         istat = 0
         return
      endif
      if ( (cv_frere(inode).eq.0) ) then
         doit=.FALSE.
         istat = 0
         return
      endif
      npiv_son2 = max(npiv/2,1)
      if( .not. MUMPS_ISTYPE2BYSIZE(nfront,npiv_son2) )then
         doit=.FALSE.
         istat = 0
         return
      endif
      ncb = nfront - npiv
      kmax = MUMPS_REG_GETKMAX(cv_keep8(21),ncb)
      if (cv_keep(50) ==  0) then
         keep48_loc=0
      else
         keep48_loc=3
      endif
      if (cv_keep(48).EQ.5) keep48_loc = 5
      if(npropmap .gt. cv_keep(83)) then
         nslaves_max   = MUMPS_BLOC2_GET_NSLAVESMAX(
     &        cv_slavef, keep48_loc, cv_keep8(21),
     &        cv_keep(50), nfront, 
     &        max(max(ncb,nfront-cv_keep(420)),1),
     &        cv_keep(375), cv_keep(119))
         nslaves_estim = min(npropmap-1,nslaves_max)
         nslaves_estim = max(nslaves_estim,1)
      else
         nslaves_max   = MUMPS_BLOC2_GET_NSLAVESMAX(
     &        cv_slavef, keep48_loc, cv_keep8(21),
     &        cv_keep(50), nfront, ncb, cv_keep(375),
     &        cv_keep(119) )
         nslaves_estim = MUMPS_BLOC2_GET_NSLAVESMIN(
     &     cv_slavef, keep48_loc,cv_keep8(21),
     &     cv_keep(50), nfront, ncb,
     &     cv_keep(375), cv_keep(119) )
         nslaves_estim = max(nslaves_estim,1)
         nslaves_estim = min(nslaves_estim,nslaves_max)
      endif
      if (cv_keep(50).eq.0) then
         wk_master = (dble(2)/dble(3))*
     &        dble(npiv)*dble(npiv)*dble(npiv)+
     &        dble(npiv)*dble(npiv)*dble(nfront-npiv)
      else
         wk_master = dble(npiv)*dble(npiv)*dble(npiv)/dble(3)
      end if
      strat = cv_keep(62)
      doit = .TRUE.
      k1 = cv_keep(82)
      k3 = cv_keep(82)
      do kk=1,cv_keep(82)-1
         npiv2   = npiv/kk
         nfront2 = nfront-npiv+npiv2
         if (npiv2 .le. max(6*cv_keep(6),0).or.
     &        (nfront2.le.cv_keep(9)) ) then
            k1 = max(1,kk-1)
            exit
         endif
         wk_master2 = wk_master / dble(kk)
         if (cv_keep(50).eq.0) then
            wk_slave2  = ( dble(npiv2)*dble(nfront2-npiv2) *
     &           dble(2*nfront2-npiv2) ) / dble(nslaves_estim)
         else
            wk_slave2 =
     &           ( dble(npiv2)*dble(nfront2-npiv2)*dble(nfront2) )
     &           /   dble(nslaves_estim)
         endif
         if(wk_master2.le.
     &        (1.0D0 +dble(kk*strat)/dble(100))*wk_slave2) then
            k1 = kk
            exit
         endif
      enddo
      do kk=1,cv_keep(82)-1
         npiv2   = npiv/kk
         nfront2 = nfront
         if (npiv2 .le. max(6*cv_keep(6),0)) then
            k3 = max(1,kk-1)
            exit
         endif
         wk_master2 = wk_master / dble(kk)
         if (cv_keep(50).eq.0) then
            wk_slave2  = ( dble(npiv2)*dble(nfront2-npiv2) *
     &           dble(2*nfront2-npiv2) ) / dble(nslaves_estim)
         else
            wk_slave2 =
     &           ( dble(npiv2)*dble(nfront2-npiv2)*dble(nfront2) )
     &           /   dble(nslaves_estim)
         endif
         if(wk_master2.le.wk_slave2) then
            k3 = kk
            exit
         endif
      enddo
      k3 = min(npiv,k3)
      k1 = min(npiv,k1)
      IF (cv_keep(79).GE.1) THEN
          k1=min(k1, npropmap-1) 
          k3=min(k3, npropmap-1) 
      ENDIF
      if(k3 .lt. k1) then
         k3 = k1
      endif
      istat=0
      return
      end subroutine MUMPS_GET_SPLIT_INKPART
      subroutine MUMPS_GET_MEMSPLIT_INKPART(inode,
     &     doit,npiv,nfront,npropmap,k2,istat)
      implicit none
      integer,intent(in)  :: inode
      logical,intent(out) :: doit
      integer,intent(in)  :: npiv,nfront,npropmap
      integer,intent(out) :: istat
      integer,intent(out) :: k2
      integer npiv2,npiv_son2
      integer kk
      DOUBLE PRECISION mem_master, mem_slave
      doit=.FALSE.
      k2=1
      istat=-1
      doit=.TRUE.
      if (cv_nodetype(inode) .gt. 0) then
         doit=.FALSE.
         istat = 0
         return
      endif
      if (cv_frere(inode).eq.0
     &   ) then
         doit=.FALSE.
         istat = 0
         return
      endif
      if ((nfront-npiv).lt.npropmap.OR.
     &    (npropmap.le.0) ) then
        doit=.FALSE.
        istat = 0
        return
      endif
      npiv_son2 = max(npiv/2,1)
      if( .not. MUMPS_ISTYPE2BYSIZE(nfront,npiv_son2) )then
        doit=.FALSE.
        istat = 0
        return
      endif
      doit = .TRUE.
      k2 = min(cv_keep(82),npropmap-1)   
      do kk=1,min(cv_keep(82)-1, npropmap-1)
         npiv2 = npiv/kk
         if(npiv2 .eq. 0) then
            k2 = max(1,kk-1)
            exit
         endif
         mem_slave   =  dble(nfront-npiv)*dble(nfront)/
     &                      dble(npropmap-kk+1)
            mem_master  = dble(npiv2)*dble(nfront)
         if(mem_master.le.
     &         (1.0D0 +dble(cv_keep(62))/dble(100))*mem_slave) then
            k2 = kk
            exit
         endif
      enddo
      k2 = max(k2, 1)
      k2 = min (k2, npiv)
      istat=0
      return
      end subroutine MUMPS_GET_MEMSPLIT_INKPART
      subroutine MUMPS_SPLITNODE_UPDATE(inode,nfront,npiv,k,
     &     lnpivsplit, npivsplit,
     &     istat)
      implicit none
      integer, intent(in)::nfront,npiv
      integer, intent(in):: k
      integer, intent(in)::lnpivsplit
      integer, intent(in)::npivsplit(lnpivsplit)
      integer, intent(in):: inode
      integer, intent(out)::istat
      integer lev,npiv_father,
     &     npiv_son,nfrontk,npivk,next_father
      DOUBLE PRECISION:: ncostm,ncostw,ncostm_ison,ncostw_ison,
     &                   ncostm_ifather,ncostw_ifather
      integer::ison,ifather
      character (len=48):: subname
      istat=-1
      subname='SPLITNODE_UPDATE'
      npiv_son = npivsplit(1)
      ison = inode
      next_father = -frere(ison)
      ncostw=cv_ncostw(inode)
      ncostm=cv_ncostm(inode)
      nfrontk = nfront
      npivk = npiv
      call MUMPS_CALCNODECOSTS(npiv_son,nfrontk,
     &     ncostw_ison,ncostm_ison)
      cv_ncostw(ison)=ncostw_ison
      cv_ncostm(ison)=ncostm_ison
      if(associated(cv_tcostw)) cv_tcostw(ison) = cv_tcostw(inode)
     &     -ncostw +cv_ncostw(ison)
      if(associated(cv_tcostm)) cv_tcostm(ison) = cv_tcostm(inode)
     &     -ncostm +cv_ncostm(ison)
      do lev = 1, k-1
         ifather = next_father
         next_father = -frere(ifather)
         npiv_son=   abs(npivsplit(lev))
         npiv_father=abs(npivsplit(lev+1))
         call MUMPS_CALCNODECOSTS(npiv_father,nfrontk-npiv_son,
     &        ncostw_ifather,ncostm_ifather)
         cv_ncostw(ifather)=ncostw_ifather
         cv_ncostm(ifather)=ncostm_ifather
         if(associated(cv_tcostw))
     &        cv_tcostw(ifather) = cv_tcostw(ison)+cv_ncostw(ifather)
         if(associated(cv_tcostm))
     &        cv_tcostm(ifather) = cv_tcostm(ison)+cv_ncostm(ifather)
         cv_total_split=cv_total_split+1
         if(lev .gt. 1) then
            call MUMPS_PROPMAP4SPLIT(inode,ison,ierr)
            if(ierr.ne.0) then
               if(cv_lp.gt.0)
     &              write(cv_lp,*)'PROPMAP4SPLIT error in ',subname
               istat = ierr
               return
            endif
         endif
         nfrontk = nfrontk-npiv_son
         npivk = npivk - npiv_son
         ison = ifather
      enddo
      if (npivk .ne. npiv_father) then
        write(*,*) "Error 1 in MUMPS_SPLITNODE_UPDATE"
        call MUMPS_ABORT()
      endif
      call MUMPS_PROPMAP4SPLIT(inode,ifather,ierr)
      if(ierr.ne.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'PROPMAP4SPLIT error in ',subname
         istat = ierr
         return
      endif
      cv_ncostw(inode) = ncostw
      cv_ncostm(inode) = ncostm
      istat = 0
      return
      end subroutine MUMPS_SPLITNODE_UPDATE
      function  MUMPS_IS_NODE_OF_TYPE2 (inode)
      implicit none
      integer, intent(in) :: inode
      logical             :: MUMPS_IS_NODE_OF_TYPE2
      if (
     &       (cv_nodetype(inode) .EQ.  2          ) .OR.
     &       (cv_nodetype(inode) .EQ.  tsplit_beg ) .OR.
     &       (cv_nodetype(inode) .EQ.  tsplit_mid ) .OR.
     &       (cv_nodetype(inode) .EQ. -tsplit_mid ) .OR.
     &       (cv_nodetype(inode) .EQ.  tsplit_last) .OR.
     &       (cv_nodetype(inode) .EQ. -tsplit_last)
     &   ) then
          MUMPS_IS_NODE_OF_TYPE2 = .TRUE.
        else
          MUMPS_IS_NODE_OF_TYPE2 = .FALSE.
        endif
      return
      end function MUMPS_IS_NODE_OF_TYPE2
      subroutine MUMPS_ENCODE_PROCNODE(istat)
      implicit none
      integer, intent(out)::istat
      integer i,in,inode
      character (len=48):: subname
      integer, external :: MUMPS_ENCODE_TPN_IPROC
      istat=-1
      subname='ENCODE_PROCNODE'
      do i=1,cv_nbsa
         inode=cv_ssarbr(i)
         cv_nodetype(inode)=0
         in=cv_fils(inode)
         do while (in>0)
            in=cv_fils(in)
         end do
         in=-in
         do while(in.gt.0)
            call MUMPS_TYPEINSSARBR(in)
            in=cv_frere(in)
         enddo
      enddo
      do i=1,cv_n
         if (cv_frere(i).lt.cv_n+1) then
            if(cv_nodetype(i).eq.cv_invalid) then
               if(cv_lp.gt.0)
     &         write(cv_lp,*)'Error in ',subname
               return
            endif
            if (i.eq.cv_keep(38)) then
               if (cv_nodetype(i).ne.3) then
                 cv_nodetype(i)=3
               endif
            endif
            cv_procnode(i)=MUMPS_ENCODE_TPN_IPROC( cv_nodetype(i),
     &           cv_procnode(i)-1, cv_keep(199))
            in=cv_fils(i)
            do while (in>0)
               cv_procnode(in)=cv_procnode(i)
               in=cv_fils(in)
            end do
         end if
      end do
      istat = 0
      return
      end subroutine MUMPS_ENCODE_PROCNODE
      subroutine MUMPS_FATHSON_REPLACE(ifather,istat)
      implicit none
      integer,intent(in)::ifather
      integer,intent(out)::istat
      integer in,son,oldl0end
      logical father_has_sons
      character (len=48):: subname
      istat=-1
      subname='FATHSON_REPLACE'
      father_has_sons=.TRUE.
      in=ifather
      do while (in.gt.0)
         in=cv_fils(in)
      end do
      if(in.eq.0) then
         cv_nodelayer(ifather)=1
         cv_keep(262)=cv_keep(262)+1
         father_has_sons=.FALSE.
      end if
      if(cv_layerl0_end-cv_layerl0_start.gt.0) then
         cv_layerl0_start= cv_layerl0_start+1
      elseif(father_has_sons) then
         cv_layerl0_start= cv_layerl0_start+1
      else
         istat=1
         cv_nodelayer(ifather)=0
         return
      endif
      cv_nbsa=cv_nbsa-1
      oldl0end = cv_layerl0_end
      if (father_has_sons) then
         son=-in
         son=-in
 10      continue
         cv_layerl0_end=cv_layerl0_end+1
         if (cv_tcostw(son).GT.mincostw)
     &           layerL0_endforarrangeL0 = layerL0_endforarrangeL0+1
         cv_layerl0_array(cv_layerl0_end)=son
         cv_layerl0_sorted_costw(cv_layerl0_end)=cv_tcostw(son)
         cv_nbsa=cv_nbsa+1
         if((cv_frere(son).gt.0).and.(cv_frere(son).lt.cv_n+1)) then
            son=cv_frere(son)
            goto 10
         end if
      endif
         cv_costw_layer0=cv_costw_layer0 - cv_ncostw(ifather)
         cv_costm_layer0=cv_costm_layer0 - cv_ncostm(ifather)
         cv_costw_upper=cv_costw_upper + cv_ncostw(ifather)
         cv_costm_upper=cv_costm_upper + cv_ncostm(ifather)
         if(cv_layerl0_end.gt.oldl0end) then
           call MUMPS_SORT_MSORT(ierr,cv_layerl0_end-oldl0end,
     &            cv_layerl0_array(oldl0end+1:cv_layerl0_end),
     &     cv_layerl0_sorted_costw(oldl0end+1:cv_layerl0_end))
           if(ierr.ne.0) then
              if(cv_lp.gt.0)
     &           write(cv_lp,*) 'Error reported by MUMPS_SORT_MSORT in',
     &           subname
              istat = ierr
              return
           endif
           call MUMPS_SORT_MMERGE(
     &        cv_layerl0_start,oldl0end,oldl0end-cv_layerl0_start+1,
     &        oldl0end+1,cv_layerl0_end,cv_layerl0_end-oldl0end,
     &        cv_layerl0_array(1:cv_layerl0_end),
     &        cv_layerl0_sorted_costw(1:cv_layerl0_end),ierr)
           if(ierr.ne.0) then
              if(cv_lp.gt.0)
     &           write(cv_lp,*)
     &           'Error reported by MUMPS_SORT_MMERGE in',
     &           subname
              istat = ierr
              return
           endif
         endif
      istat=0
      return
      end subroutine MUMPS_FATHSON_REPLACE
      subroutine MUMPS_FIND_BEST_PROC(inode,map_strat,work,mem,
     &                       workload,memused,proc,istat,respect_prop)
!DEC$ NOOPTIMIZE
      implicit none
      integer, intent(in)::inode,map_strat
      DOUBLE PRECISION,intent(in)::work,mem
      DOUBLE PRECISION,dimension(:),intent(inout)::workload, memused
      integer,intent(out):: proc,istat
      logical,intent(in),OPTIONAL::respect_prop
      integer i
      logical respect_proportional
      intrinsic huge
      DOUBLE PRECISION dummy
      character (len=48):: subname
      istat=-1
      respect_proportional=.FALSE.
      if(present(respect_prop)) respect_proportional=respect_prop
      subname='FIND_BEST_PROC'
      proc=-1
      if((map_strat.ne.cv_equilib_flops).and.
     &   (map_strat.ne.cv_equilib_mem)) return
      dummy=huge(dummy) 
      do i=cv_slavef,1,-1
         if (
     &       ((.NOT.respect_proportional)
     &        .OR.
     &        (MUMPS_BIT_GET4PROC(inode,i).AND.respect_proportional))
     &      .AND.
     &       (((workload(i).lt.dummy).AND.
     &                         (map_strat.eq.cv_equilib_flops))
     &        .OR.
     &       ((memused(i).lt.dummy).AND.
     &                         (map_strat.eq.cv_equilib_mem))))then
            if((.not.cv_constr_work).or.
     &         (workload(i)+work.lt.cv_proc_maxwork(i))) then
               if((.not.cv_constr_mem).or.
     &            (memused(i)+mem.lt.cv_proc_maxmem(i))) then
                  proc=i
                  if(map_strat.eq.cv_equilib_flops) then
                     dummy=workload(i)
                  elseif(map_strat.eq.cv_equilib_mem) then
                     dummy=memused(i)
                  endif
               end if
            end if
         end if
      end do
      if (proc.ne.-1) then
         workload(proc)=workload(proc)+work
         memused(proc)=memused(proc)+mem
         istat=0
      end if
      return
      end subroutine MUMPS_FIND_BEST_PROC
      subroutine MUMPS_FIND_THISLAYER(nmb,
     &                                   thislayer,nmb_thislayer,istat)
      implicit none
      integer, intent(in)::nmb
      integer,intent(out) :: thislayer(:) 
      integer,intent(out) :: nmb_thislayer,istat
      integer i
      character (len=48):: subname
      istat=-1
      subname='FIND_THISLAYER'
      thislayer=0
      nmb_thislayer=0
      if((nmb.lt.0).or.(nmb.gt.cv_maxlayer)) return
      do i=1,cv_n
         if(cv_nodelayer(i).eq.nmb) then
            nmb_thislayer=nmb_thislayer+1
            if(nmb_thislayer.gt.cv_maxnodenmb) then
               if(cv_lp.gt.0)
     &         write(cv_lp,*)'Problem with nmb_thislayer in ',subname
               return
            endif
            thislayer(nmb_thislayer)=i
         end if
      end do
      istat=0
      return
      end subroutine MUMPS_FIND_THISLAYER
      subroutine MUMPS_HIGHER_LAYER(startlayer,thislayer,
     &                 nmb_thislayer,cont,istat)
      implicit none
      integer,intent(in)::startlayer,nmb_thislayer
      integer,intent(in)::thislayer(:) 
      logical,intent(inout)::cont
      integer,intent(out)::istat 
      integer :: visited
      integer il,i,current,in,ifather
      logical father_valid,upper_layer_exists
      character (len=48):: subname
      istat=-1
      subname='HIGHER_LAYER'
      if(.NOT.cont) return
      if(startlayer.lt.1) return
      current=startlayer-1
      visited = -current-1
         upper_layer_exists=.FALSE.
         if (current.eq.0) then
          do i=1,cv_n
            if (cv_nodelayer(i).ne.current) then
               if(cv_nodelayer(i).eq.1) then
                  upper_layer_exists=.TRUE.
                  exit
               endif
            endif
          enddo
         endif
         do il=1,nmb_thislayer
            i = thislayer(il)
            in=i
            if (cv_nodetype(in).eq.tsplit_beg) then
              do while (cv_frere(in).lt.0)
                ifather = -cv_frere(in)
                if  (abs(cv_nodetype(ifather)).eq.tsplit_mid) then
                  in  = ifather 
                  cv_nodelayer (in) = -visited-1
                  cycle
                else if  (abs(cv_nodetype(ifather)).eq.tsplit_last) then
                  in = ifather
                  cv_nodelayer (in) = current
                  exit
               else
                 write(6,*) ' Internal error 1 in MUMPS_HIGHER_LAYER'
                 call MUMPS_ABORT()
                endif
              end do
            endif
         enddo
         do il=1,nmb_thislayer
            i = thislayer(il)
            if (cv_nodelayer(i).lt.current) cycle
            in=i
            if (cv_nodetype(in).eq.tsplit_beg) then
              cv_nodelayer (in) = visited
              do while (cv_frere(in).lt.0)
                ifather = -cv_frere(in)
                if  (abs(cv_nodetype(ifather)).eq.tsplit_mid) then
                  in  = ifather 
                  cv_nodelayer (in) = -visited-1
                  cycle
                else if (abs(cv_nodetype(ifather)).eq.tsplit_last) then
                  in = ifather
                  exit
               else
                 write(6,*) ' Internal error 1 in MUMPS_HIGHER_LAYER',
     &           cv_nodetype(ifather)
                 call MUMPS_ABORT()
                endif
              end do
            endif
            if(cv_frere(in).eq.0) cycle
            cv_nodelayer (in) = visited
               father_valid=.TRUE.
            do while(cv_frere(in).gt.0)
               if (cv_nodelayer(cv_frere(in)).gt.current) then
                  father_valid=.FALSE.
                  in = cv_frere(in)
                  cycle
               endif
               if (cv_nodelayer(cv_frere(in)).eq.visited) exit
               in=cv_frere(in) 
               if (cv_nodelayer(in).eq.current)  then
                    cv_nodelayer(in) = visited
               endif
            end do
            if (.not.father_valid .or. cv_frere(in).gt.0) then 
                    cycle
            endif
            ifather=-cv_frere(in)
            if(cv_nodelayer(ifather).eq.current+1) then
               cycle
            endif
            in=ifather
            do while (cv_fils(in).gt.0)
               in=cv_fils(in)
            end do
            in=-cv_fils(in)
            if(cv_nodelayer(in).gt.current) then 
               father_valid=.FALSE.
            else
               father_valid=.TRUE.
               do while(cv_frere(in).gt.0)
                in=cv_frere(in)
                if(cv_nodelayer(in).gt.current) then 
                       father_valid=.FALSE.
                       exit
                endif
                if(cv_nodelayer(in).eq.visited) then
                    exit  
                endif
               end do
            endif
            if(father_valid) then
               cv_nodelayer(ifather)=current+1
               upper_layer_exists=.TRUE.
            end if
         end do
         if (upper_layer_exists) then 
             current=current+1
             cv_maxlayer=current
             cont=.TRUE.
         else
             cv_maxlayer=current
             cont=.FALSE.
         endif
         do il=1,nmb_thislayer
            i = thislayer(il)
           if (cv_nodelayer(i).eq.visited) cv_nodelayer(i) = -visited-1
         enddo
      istat=0
      return
      end subroutine MUMPS_HIGHER_LAYER
      subroutine MUMPS_INITPART1(n,slavef,
     &                    frere,fils,nfsiz,ne,keep,KEEP8,icntl,info,
     &                    procnode,ssarbr,peak,istat
     &                      , SIZEOFBLOCKS, LSIZEOFBLOCKS
     &     )
      implicit none
      integer, intent(in)::n,slavef
      integer, intent(in), TARGET:: frere(n),fils(n),nfsiz(n),ne(n),
     &  keep(500),icntl(60),info(80),
     &  procnode(n),ssarbr(n)
      INTEGER(8), intent(in), TARGET:: KEEP8(150)
      integer,intent(out)::istat
      integer, intent(in)         :: LSIZEOFBLOCKS
      integer, intent(in), TARGET :: SIZEOFBLOCKS(LSIZEOFBLOCKS)
      integer i,allocok,rest
      DOUBLE PRECISION peak
      character (len=48):: subname
      intrinsic bit_size,min,max
      istat=-1
      nullify(cv_frere,cv_fils,cv_nfsiz,cv_ne,cv_keep,cv_keep8,
     &        cv_icntl,cv_info,cv_procnode,cv_ssarbr)
      nullify(cv_ncostw,cv_tcostw,cv_ncostm,cv_tcostm,
     &        cv_nodelayer,cv_nodetype,cv_depth,
     &        cv_layerworkload,cv_layermemused,cv_prop_map)
       nullify(cv_SIZEOFBLOCKS)
       cv_SIZEOFBLOCKS => SIZEOFBLOCKS
      subname='INITPART1'
      cv_n=n
      cv_slavef=slavef
      cv_keep=>keep
      cv_keep8=>KEEP8
      if(cv_keep(82) .lt. 0) then
         write(cv_lp,*)
     &        'Warning in mumps_static_mapping : splitting is set off'
         cv_keep(82) = 0
      endif
      if(cv_keep(83) .lt. 0) then
         write(cv_lp,*)
     &        'warning in mumps_static_mapping : keep(83) reset to 0'
         cv_keep(83) = 0
      endif
      if(slavef.gt.1) then
         cv_mixed_strat_bound = max(cv_keep(78),1)
         cv_maxdepth = slavef
      else
         cv_maxdepth = 0
         cv_mixed_strat_bound=0
      endif
      cv_bitsize_of_int = bit_size(n)
      if(cv_bitsize_of_int.le.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Problem with bit size in ',subname
         return
      endif
      rest = mod(cv_slavef,cv_bitsize_of_int)
      if (rest.eq.0) then
         cv_size_ind_proc = cv_slavef / cv_bitsize_of_int
      else
         cv_size_ind_proc = cv_slavef / cv_bitsize_of_int + 1
      endif
      allocate(cv_ncostw(n),cv_tcostw(n),cv_ncostm(n),cv_tcostm(n),
     &       cv_nodelayer(n),cv_nodetype(n),cv_depth(n),
     &       cv_layerworkload(slavef),cv_layermemused(slavef),
     &       cv_prop_map(n),STAT=allocok)
      if (allocok.gt.0) then
         cv_info(1) = cv_error_memalloc
         cv_info(2) = 8*n+2*cv_slavef
         istat = cv_error_memalloc
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'memory allocation error in ',subname
         return
      end if
      if(cv_keep(82) .eq. 0) then
            if(cv_lp.gt.0)
     &           write(cv_lp,*)' No splitting during static mapping '
      endif
      cv_frere=>frere
      cv_fils=>fils
      cv_nfsiz=>nfsiz
      cv_ne=>ne
      cv_icntl=>icntl
      cv_info=>info
      cv_procnode=>procnode
      cv_ssarbr=>ssarbr
      cv_ssarbr=0
      cv_nodetype=cv_invalid
      cv_nsteps=keep(28)
      if((keep(28).gt.n).OR.(keep(28).lt.0)) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'problem with nsteps in ',subname
         return
      end if
      cv_costw_upper=0.0D0
      cv_costm_upper=0.0D0
      cv_costw_layer0=0.0D0
      cv_costm_layer0=0.0D0
      cv_costw_total=0.0D0
      cv_costm_total=0.0D0
      cv_nodelayer=n+2
      cv_depth=cv_invalid
      cv_l0wthresh=0.0D0 
      cv_splitthresh=0.45D0
      cv_relax=dble(1) + dble(max(0,keep(68)))/dble(100)
      cv_maxlayer=0
      cv_maxnsteps= cv_nsteps+1
      cv_layerworkload=dble(0)
      cv_layermemused=dble(0)
      cv_total_amalg=0
      cv_total_split=0
      cv_last_splitting%new_ison=cv_invalid
      cv_last_splitting%new_ifather=cv_invalid
      cv_last_splitting%old_keep2=cv_invalid
      cv_last_splitting%ncostw_oldinode=cv_d_invalid
      cv_last_splitting%ncostm_oldinode=cv_d_invalid
      cv_last_splitting%tcostw_oldinode=cv_d_invalid
      cv_last_splitting%tcostm_oldinode=cv_d_invalid
      do i=1,cv_n
         nullify(cv_prop_map(i)%ind_proc)
      end do
      istat=0
      return
      end subroutine MUMPS_INITPART1
      subroutine MUMPS_INITPART2(istat)
      implicit none
      integer,intent(out)::istat
      integer i,allocok,inode,in,inoderoot,ierr,maxcut
      character (len=48):: subname
      istat=-1
      subname='INITPART2'
      if(associated(cv_layerl0_array))deallocate(cv_layerl0_array)
      if(associated(cv_layerl0_sorted_costw))
     &           deallocate(cv_layerl0_sorted_costw)
      deallocate(cv_depth,cv_tcostw,cv_tcostm,STAT=ierr)
      if(ierr.ne.0) then
         if(cv_lp.gt.0)
     &   write(cv_lp,*)'Memory deallocation error in ',subname
         istat = cv_error_memdeloc
         return
      end if
      if(cv_maxnsteps.lt.1) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'problem with maxnsteps in ',subname
         return
      end if
      cv_maxnodenmb=cv_maxnsteps
      do i=1,cv_nbsa
         inode=cv_ssarbr(i)
         inoderoot=inode
 300     continue
         in = inode
         do while (in.ne.0)
            inode = in
            do while (in.gt.0)
            in = cv_fils(in)
            end do
            if (in.lt.0) in=-in
         end do
 100     continue
         if (inode.ne.inoderoot) then
            cv_maxnodenmb=cv_maxnodenmb-1
            in = cv_frere(inode)
            inode = abs(in)
            if (in.lt.0) then
               go to 100
            else
               go to 300
            end if
         end if
      end do
      if(cv_keep(82) .gt. 0) then
         maxcut = min((cv_keep(82)-1)*cv_maxnodenmb,cv_n)
         cv_maxnsteps = min(cv_maxnsteps+maxcut,cv_n)
         cv_maxnodenmb = min(cv_maxnodenmb+maxcut,cv_n)
      endif
      nullify(cv_layer_p2node)
      if(cv_maxnodenmb.lt.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'problem with maxnodenmb in ',subname
         return
      elseif(cv_maxnodenmb.lt.1) then
         cv_maxnodenmb = 1
      end if
      allocate(cv_layer_p2node(cv_maxnodenmb),STAT=allocok)
      if (allocok.gt.0) then
         cv_info(1) = cv_error_memalloc
         cv_info(2) = cv_maxnodenmb
         istat = cv_error_memalloc
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'memory allocation error in ',subname
         return
      end if
      do i=1,cv_maxnodenmb
         nullify(cv_layer_p2node(i)%t2_nodenumbers,
     &           cv_layer_p2node(i)%t2_cand,
     &           cv_layer_p2node(i)%t2_candcostw,
     &           cv_layer_p2node(i)%t2_candcostm)
         cv_layer_p2node(i)%nmb_t2s=0
      enddo
      istat = 0
      end subroutine MUMPS_INITPART2
      function MUMPS_ISTYPE2BYSIZE(nfront,npiv)
      implicit none
      logical::MUMPS_ISTYPE2BYSIZE
      integer,intent(in)::nfront,npiv
      MUMPS_ISTYPE2BYSIZE=.FALSE.
      if( (nfront - npiv > cv_keep(9))
     &     .and. ((npiv > cv_keep(4)).or.(.TRUE.))
     &     .and. (cv_icntl(59).eq.0) ) MUMPS_ISTYPE2BYSIZE=.TRUE.
      return
      end function MUMPS_ISTYPE2BYSIZE
      subroutine MUMPS_LAYERL0(istat)
      implicit none
      integer,intent(out)::istat
      integer i,ierr,inode
      logical accepted
      integer,parameter::map_strat=cv_equilib_flops
      character (len=48):: err_rep,subname
      logical use_geist_ng_replace, skiparrangeL0
      INTEGER MINSIZE_L0
      INTEGER CURRENT_SIZE_L0
      istat=-1
      subname='LAYERL0'
      accepted=.FALSE.
      IF (cv_keep(72).EQ.2) THEN
       MINSIZE_L0 = 6*cv_slavef
      ELSE
       IF (cv_keep(66).NE.0) THEN
         IF (cv_keep(66).EQ.1) THEN
           MINSIZE_L0 = 3*cv_slavef
         ELSE
           MINSIZE_L0 = 2*cv_slavef
         ENDIF
       ELSE
         MINSIZE_L0 = 3*cv_slavef
       ENDIF
      ENDIF
 55   continue
      skiparrangeL0 = .false.
      do while(.not.accepted)
         IF (cv_keep(66).EQ.2) THEN
            CURRENT_SIZE_L0 = layerL0_endforarrangeL0
         ELSE
            CURRENT_SIZE_L0 = layerL0_endforarrangeL0
         ENDIF
         IF ( (    (CURRENT_SIZE_L0.LT.MINSIZE_L0)
     &             .OR. skiparrangeL0
     &        )
     &        .AND.
     &           (cv_layerl0_end.LT.cv_maxnsteps/2) ) THEN
          accepted = .false.  
         ELSE
          err_rep='ARRANGEL0'
          call MUMPS_ARRANGEL0(map_strat, layerL0_endforarrangeL0,
     &                           cv_layerworkload,cv_layermemused,
     &                           cv_procnode,ierr)
          if(ierr.ne.0) then
            if(cv_lp.gt.0)
     &      write(cv_lp,*)'Error reported by ',err_rep,' in ',subname
            istat = ierr
            return
          end if
          err_rep='ACCEPT_L0'
          call MUMPS_ACCEPT_L0(map_strat,
     &                           cv_layerworkload,cv_layermemused,
     &                           accepted,ierr)
          if(ierr.ne.0) then
            if(cv_lp.gt.0)
     &      write(cv_lp,*)'Error reported by ',err_rep,' in ',subname
            istat = ierr
            return
          end if
         ENDIF
         IF (cv_keep(66).EQ.0) THEN
           IF (cv_slavef.GT.16) 
     &          skiparrangeL0 = .NOT.skiparrangeL0   
         ENDIF
         if (accepted.OR.(cv_costw_total.le.0.0D0)) then
            exit
         elseif(((cv_costw_layer0/cv_costw_total).gt.cv_l0wthresh) .AND.
     &           (.TRUE.))then 
            err_rep='MAX_TCOST_L0'
            inode = cv_layerl0_array(cv_layerl0_start)
            use_geist_ng_replace = .TRUE.
            if(use_geist_ng_replace) then
               err_rep='FATHSON_REPLACE'
               call MUMPS_FATHSON_REPLACE(inode,ierr)
               if(ierr.eq.1) then
                  accepted=.TRUE.
               elseif(ierr.ne.0) then
                  if(cv_lp.gt.0)
     &            write(cv_lp,*)
     &            'Error rep. by ',err_rep,' in ',subname
                  istat = ierr
                  return
               endif
            endif
         else
            accepted=.TRUE.
         end if
      end do
      accepted=.TRUE.
      if (accepted) then
      else
         goto 55
      endif
      err_rep='LIST2LAYER'
      call MUMPS_LIST2LAYER(ierr)
      if(ierr.ne.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error reported by ',err_rep,' in ',subname
         istat = ierr
         return
      end if
      err_rep='MAKE_PROPMAP'
      call MUMPS_MAKE_PROPMAP(ierr)
      if(ierr.ne.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error reported by ',err_rep,' in ',subname
         istat = ierr
         return
      end if
      if ( cv_keep(75).EQ.1 ) then
         call MUMPS_ARRANGEL0(map_strat, cv_layerl0_end,
     &                     cv_layerworkload,cv_layermemused,
     &                     cv_procnode,ierr, respect_prop=.TRUE.)
         if(ierr.ne.0) then
            if(cv_lp.gt.0)
     &      write(cv_lp,*)'Error reported by ',err_rep,' in ',subname
            istat = ierr
            return
         end if
      else if (layerL0_endforarrangeL0.LT.cv_layerl0_end) THEN
         call MUMPS_ARRANGEL0(map_strat, cv_layerl0_end,
     &                     cv_layerworkload,cv_layermemused,
     &                     cv_procnode,ierr)
      endif
      call MUMPS_MAPSUBTREE(cv_procnode)
      do i=1,cv_slavef
         cv_proc_workload(i)=cv_layerworkload(i)
         cv_proc_memused(i)=cv_layermemused(i)
      end do
      istat=0
      return
      end subroutine MUMPS_LAYERL0
      subroutine MUMPS_LIST2LAYER(istat)
      implicit none
      integer, intent(out)::istat
      character (len=48):: subname
      integer i,inode
      istat=-1
      subname='LIST2LAYER'
      cv_dist_L0_mixed_strat_bound=0
      cv_nbsa=0
      do i=cv_layerl0_start,cv_layerl0_end
         inode=cv_layerl0_array(i)
         if(inode.gt.0) then
         cv_dist_L0_mixed_strat_bound=max(cv_dist_L0_mixed_strat_bound
     &        ,max(cv_depth(inode)-cv_mixed_strat_bound,0))
         cv_nodelayer(inode)=0
         cv_nbsa=cv_nbsa+1
         cv_ssarbr(cv_nbsa)=inode
         endif
      enddo
      istat=0
      return
      end subroutine MUMPS_LIST2LAYER
      subroutine MUMPS_MAKE_PROPMAP(istat)
      implicit none
      integer,intent(out)::istat
      integer i,pctr,pctr2,ierr
      character (len=48):: subname
      INTEGER, ALLOCATABLE, DIMENSION(:) :: procindex
      INTEGER :: allocok
      subname = "MUMPS_MAKE_PROPMAP"
      istat = -1
      ALLOCATE(procindex(cv_size_ind_proc),stat=allocok)
      IF (allocok > 0) THEN
        cv_info(1) = cv_error_memalloc
        cv_info(2) = cv_maxnodenmb
        istat = cv_error_memalloc
        if(cv_lp.gt.0)
     &       write(cv_lp,*) 'Memory allocation error in ',subname
        return
      ENDIF
      pctr=cv_n
      pctr2=cv_mixed_strat_bound
      do i=1,cv_slavef
         call MUMPS_BIT_SET(procindex,i,ierr)
         if(ierr.ne.0) then
            if(cv_lp.gt.0)write(cv_lp,*)
     &           'MUMPS_BIT_SET signalled error to',subname
            istat = ierr
            GOTO 999
         end if
      end do
      do i=1,cv_n
         if(cv_frere(i).eq.0) then
            if(.NOT.associated(cv_prop_map(i)%ind_proc)) then
               call MUMPS_PROPMAP_INIT(i,ierr)
               if(ierr.ne.0) then
                  if(cv_lp.gt.0)
     &                 write(cv_lp,*)'PROPMAP_INIT signalled error to'
     &                 ,subname
                  istat = ierr
                  GOTO 999
               end if
            endif
            cv_prop_map(i)%ind_proc = procindex
            call MUMPS_PROPMAP(i,pctr,ierr)
            if(ierr.ne.0) then
            if(cv_lp.gt.0)write(cv_lp,*)
     &           'PROPMAP signalled error to',subname
               istat = ierr
               GOTO 999
            endif
            if((cv_keep(24).eq.16).OR.(cv_keep(24).eq.18)) then
               call MUMPS_MOD_PROPMAP(i,pctr2,ierr)
               if(ierr.ne.0) then
                  if(cv_lp.gt.0)write(cv_lp,*)
     &           'MOD_PROPMAP signalled error to',subname
                  istat = ierr
                  GOTO 999
               endif
            endif
         endif
      end do
      istat = 0
 999  CONTINUE
      DEALLOCATE(procindex)
      return
      end subroutine MUMPS_MAKE_PROPMAP
      subroutine MUMPS_MAP_LAYER(layernmb,thislayer,
     &   nmb_thislayer,map_strat,istat)
      implicit none
      integer, intent(in)::layernmb,thislayer(:), 
     &                     nmb_thislayer,map_strat
      integer,intent(out)::istat
      integer i,inode,j,k,ierr,nmb,aux_int,nmb_cand_needed
      DOUBLE PRECISION aux_flop,aux_mem
      INTEGER, ALLOCATABLE, DIMENSION(:) :: candid, sorted_nmb
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: 
     &  sorted_costw, sorted_costm, old_workload, old_memused
      character (len=48):: err_rep,subname
      logical use_propmap
      istat=-1
      subname='MAP_LAYER'
      if((cv_keep(24).eq.8).OR.(cv_keep(24).eq.10)
     &   .OR.(cv_keep(24).eq.12).OR.(cv_keep(24).eq.14)
     &   .OR.(cv_keep(24).eq.16).OR.(cv_keep(24).eq.18)) then
         use_propmap=.TRUE.
      else
         use_propmap=.FALSE.
      endif
      if((layernmb.lt.0).or.(layernmb.gt.cv_maxlayer)) return
      if((map_strat.ne.cv_equilib_flops).and.
     &   (map_strat.ne.cv_equilib_mem)) return
      ALLOCATE(candid(cv_slavef), sorted_nmb(2*nmb_thislayer),
     & sorted_costw(2*nmb_thislayer), sorted_costm(2*nmb_thislayer),
     & old_workload(cv_slavef), old_memused(cv_slavef), stat=allocok)
      if (allocok.gt.0) then
         cv_info(1) = cv_error_memalloc
         cv_info(2) = 7*nmb_thislayer+2*cv_slavef
         istat = cv_error_memalloc
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'memory allocation error in ',subname
         goto 999
      end if
      do i=1,nmb_thislayer
         inode=thislayer(i)
         if (cv_nodetype(inode).eq.3) then
            cv_procnode(inode)=1
            exit
         end if
      end do
      do i=1,cv_slavef
         old_workload(i)=cv_layerworkload(i)
         old_memused(i)=cv_layermemused(i)
      enddo
      nmb=0
      do i=1,nmb_thislayer
         inode=thislayer(i)
         if(cv_nodetype(inode).eq.1) then
            nmb=nmb+1
            sorted_nmb(nmb)=inode
            sorted_costw(nmb)=cv_ncostw(inode)
            sorted_costm(nmb)=cv_ncostm(inode)
         else if(MUMPS_IS_NODE_OF_TYPE2(inode)) then
            nmb=nmb+1
            do j=1,cv_layer_p2node(layernmb)%nmb_t2s
               if(cv_layer_p2node(layernmb)%t2_nodenumbers(j).ne.inode)
     &            then
                  cycle
               else
                  sorted_costw(nmb)=
     &                 cv_layer_p2node(layernmb)%t2_candcostw(j)
                  sorted_costm(nmb)=
     &                 cv_layer_p2node(layernmb)%t2_candcostm(j)
               endif
            enddo
            if((sorted_costw(nmb).eq.cv_d_invalid).OR.
     &           (sorted_costm(nmb).eq.cv_d_invalid)) then
               if(cv_lp.gt.0)
     &         write(cv_lp,*)'Error in ',subname
               goto 999
            end if
            if(sorted_costw(nmb).lt.cv_ncostw(inode))then
               sorted_costw(nmb)=cv_ncostw(inode)
               sorted_costm(nmb)=cv_ncostm(inode)
               sorted_nmb(nmb)=inode
            else
               sorted_nmb(nmb)=-inode
            endif
         else if(cv_nodetype(inode).eq.3) then
            cycle
         else
            if(cv_lp.gt.0)
     &      write(cv_lp,*)'Unknown node type. Error in ',subname
            goto 999
         end if
      end do
      if (map_strat.eq.cv_equilib_flops) then
         call MUMPS_SORT_MSORT(ierr,nmb,sorted_nmb(1:nmb),
     &                   sorted_costw(1:nmb),sorted_costm(1:nmb))
      elseif(map_strat.eq.cv_equilib_mem) then
         call MUMPS_SORT_MSORT(ierr,nmb,sorted_nmb(1:nmb),
     &                   sorted_costm(1:nmb),sorted_costw(1:nmb))
      endif
      if(ierr.ne.0) then
         if(cv_lp.gt.0)
     &   write(cv_lp,*)
     &   'Error reported by MUMPS_SORT_MSORT in ',subname
         istat = ierr
         goto 999
      endif
      do i=1,nmb
         aux_int=sorted_nmb(i)
         aux_flop=sorted_costw(i)
         aux_mem=sorted_costm(i)
         k=1
         if (aux_int.lt.0) then
            inode=-aux_int
            err_rep='SORTPROCS'
            if(use_propmap) then
               call MUMPS_SORTPROCS(map_strat,
     &              cv_proc_workload,cv_proc_memused,
     &              inode=inode,istat=ierr)
            else
               call MUMPS_SORTPROCS(map_strat,
     &              cv_proc_workload,cv_proc_memused,
     &              istat=ierr)
            end if
            if(ierr.ne.0) then
               if(cv_lp.gt.0)
     &         write(cv_lp,*)
     &         'Error reported by ',err_rep,' in ',subname
               istat = ierr
               goto 999
            endif
            nmb_cand_needed=cv_invalid
            do j=1,cv_layer_p2node(layernmb)%nmb_t2s
               if(cv_layer_p2node(layernmb)%t2_nodenumbers(j).ne.inode)
     &            then
                  cycle
               else
                  nmb_cand_needed=
     &                  cv_layer_p2node(layernmb)%t2_cand(j,cv_slavef+1)
                  exit
               endif
            enddo
            if(nmb_cand_needed.eq.cv_invalid) then
               if(cv_lp.gt.0)
     &         write(cv_lp,*)'Error in ',subname
               goto 999
            endif
            do while((k.le.cv_slavef).and.(nmb_cand_needed.gt.0))
               if(((.not.cv_constr_work).or.
     &            (cv_proc_workload(cv_proc_sorted(k))+aux_flop.lt.
     &             cv_proc_maxwork(cv_proc_sorted(k))))
     &           .AND.((.not.cv_constr_mem).or.
     &            (cv_proc_memused(cv_proc_sorted(k))+aux_mem.lt.
     &             cv_proc_maxmem(cv_proc_sorted(k))))
     &           .AND. 
     & (cv_layer_p2node(layernmb)%t2_cand(j,cv_proc_sorted(k)).eq.0))
     &            then
                  cv_proc_workload(cv_proc_sorted(k))=
     &                 cv_proc_workload(cv_proc_sorted(k))+aux_flop
                  cv_proc_memused(cv_proc_sorted(k))=
     &                 cv_proc_memused(cv_proc_sorted(k))+aux_mem
                  cv_layer_p2node(layernmb)%t2_cand(j,cv_proc_sorted(k))
     &                                             =inode
                  cv_layerworkload(cv_proc_sorted(k))=
     &                 cv_layerworkload(cv_proc_sorted(k))+aux_flop
                  cv_layermemused(cv_proc_sorted(k))=
     &                 cv_layermemused(cv_proc_sorted(k))+aux_mem
                  nmb_cand_needed=nmb_cand_needed-1
                  k=k+1
               else
                  k=k+1
                  if(k.gt.cv_slavef) then
                     if(cv_lp.gt.0)
     &                    write(cv_lp,*)'Error in ',subname
                     goto 999
                  endif
               end if
            end do
            if(nmb_cand_needed.gt.0) then
               if(cv_lp.gt.0)
     &         write(cv_lp,*)'Error in ',subname
               goto 999
            endif
            aux_flop=cv_ncostw(inode)
            aux_mem=cv_ncostm(inode)
            do while(k.le.cv_slavef)
               if(((.not.cv_constr_work).or.
     &              (cv_proc_workload(cv_proc_sorted(k))+aux_flop.lt.
     &              cv_proc_maxwork(cv_proc_sorted(k))))
     &              .AND.((.not.cv_constr_mem).or.
     &              (cv_proc_memused(cv_proc_sorted(k))+aux_mem.lt.
     &              cv_proc_maxmem(cv_proc_sorted(k))))
     &              .AND.       
     & (cv_layer_p2node(layernmb)%t2_cand(j,cv_proc_sorted(k)).eq.0))
     &              then
                  cv_procnode(inode)=cv_proc_sorted(k)
                  cv_proc_workload(cv_proc_sorted(k))=
     &                 cv_proc_workload(cv_proc_sorted(k))+aux_flop
                  cv_proc_memused(cv_proc_sorted(k))=
     &                 cv_proc_memused(cv_proc_sorted(k))+aux_mem
                  cv_layer_p2node(layernmb)%t2_cand(j,cv_proc_sorted(k))
     &                 =-inode
                  cv_layerworkload(cv_proc_sorted(k))=
     &                 cv_layerworkload(cv_proc_sorted(k))+aux_flop
                  cv_layermemused(cv_proc_sorted(k))=
     &                 cv_layermemused(cv_proc_sorted(k))+aux_mem
                  exit
               else
                  k=k+1
                  if(k.gt.cv_slavef) then
                     if(cv_lp.gt.0)
     &                    write(cv_lp,*)'Error in ',subname
                     goto 999
                  endif
               end if
            end do
         else
            inode=aux_int
            err_rep='SORTPROCS'
            if(use_propmap) then
               call MUMPS_SORTPROCS(map_strat,
     &                        cv_proc_workload,cv_proc_memused,
     &                        inode=inode,istat=ierr)
            else
               call MUMPS_SORTPROCS(map_strat,
     &                        cv_proc_workload,cv_proc_memused,
     &                        inode,istat=ierr)
            endif
            if(ierr.ne.0) then
               if(cv_lp.gt.0)
     &          write(cv_lp,*)
     &          'Error reported by ',err_rep,' in ',subname
               istat = ierr
               goto 999
            endif
            if (cv_nodetype(inode).eq.1) then
               do while(k.le.cv_slavef)
                  if((.not.cv_constr_work).or.
     &               (cv_proc_workload(cv_proc_sorted(k))+aux_flop.lt.
     &                cv_proc_maxwork(cv_proc_sorted(k)))
     &            .AND.((.not.cv_constr_mem).or.
     &               (cv_proc_memused(cv_proc_sorted(k))+aux_mem.lt.
     &                cv_proc_maxmem(cv_proc_sorted(k))))) then
                     cv_procnode(inode)=cv_proc_sorted(k)
                     cv_proc_workload(cv_proc_sorted(k))=
     &                    cv_proc_workload(cv_proc_sorted(k))+aux_flop
                     cv_proc_memused(cv_proc_sorted(k))=
     &                    cv_proc_memused(cv_proc_sorted(k))+aux_mem
                     cv_layerworkload(cv_proc_sorted(k))=
     &                    cv_layerworkload(cv_proc_sorted(k))+aux_flop
                     cv_layermemused(cv_proc_sorted(k))=
     &                    cv_layermemused(cv_proc_sorted(k))+aux_mem
                     exit
                  else
                     k=k+1
                     if(k.gt.cv_slavef) then
                        if(cv_lp.gt.0)
     &                  write(cv_lp,*)'Inconsist data in ',subname
                        goto 999
                     endif
                  end if
               end do
            elseif (MUMPS_IS_NODE_OF_TYPE2(inode)) then
               do j=1,cv_layer_p2node(layernmb)%nmb_t2s
                  if(cv_layer_p2node(layernmb)%t2_nodenumbers(j).ne.
     &               inode) then
                     cycle
                  else
                     exit
                  endif
               enddo
               do while(k.le.cv_slavef)
                  if(((.not.cv_constr_work).or.
     &               (cv_proc_workload(cv_proc_sorted(k))+aux_flop.lt.
     &                cv_proc_maxwork(cv_proc_sorted(k))))
     &              .AND.((.not.cv_constr_mem).or.
     &               (cv_proc_memused(cv_proc_sorted(k))+aux_mem.lt.
     &                cv_proc_maxmem(cv_proc_sorted(k))))
     &              .AND. 
     & (cv_layer_p2node(layernmb)%t2_cand(j,cv_proc_sorted(k)).eq.0))
     &               then
                     cv_procnode(inode)=cv_proc_sorted(k)
                     cv_proc_workload(cv_proc_sorted(k))=
     &                    cv_proc_workload(cv_proc_sorted(k))+aux_flop
                     cv_proc_memused(cv_proc_sorted(k))=
     &                     cv_proc_memused(cv_proc_sorted(k))+aux_mem
                cv_layer_p2node(layernmb)%t2_cand(j,cv_proc_sorted(k))
     &                                           =-inode
                     cv_layerworkload(cv_proc_sorted(k))=
     &                    cv_layerworkload(cv_proc_sorted(k))+aux_flop
                     cv_layermemused(cv_proc_sorted(k))=
     &                    cv_layermemused(cv_proc_sorted(k))+aux_mem
                     exit
                  else
                     k=k+1
                     if(k.gt.cv_slavef) then
                        if(cv_lp.gt.0)
     &                       write(cv_lp,*)'Error in ',subname
                        goto 999
                     endif
                  end if
               end do
               nmb_cand_needed=cv_invalid
               do j=1,cv_layer_p2node(layernmb)%nmb_t2s
                  if(cv_layer_p2node(layernmb)%t2_nodenumbers(j)
     &                 .ne.inode)
     &                 then
                     cycle
                  else
                     nmb_cand_needed=
     &                    cv_layer_p2node(layernmb)%
     &                    t2_cand(j,cv_slavef+1)
                     exit
                  endif
               enddo
               if(nmb_cand_needed.eq.cv_invalid) then
                  if(cv_lp.gt.0)
     &                 write(cv_lp,*)'Error in ',subname
                  goto 999
               endif
               aux_flop=
     &              cv_layer_p2node(layernmb)%t2_candcostw(j)
               aux_mem=
     &              cv_layer_p2node(layernmb)%t2_candcostm(j)
               do while((k.le.cv_slavef).and.(nmb_cand_needed.gt.0))
                  if(((.not.cv_constr_work).or.
     &                 (cv_proc_workload(cv_proc_sorted(k))+aux_flop.lt.
     &                 cv_proc_maxwork(cv_proc_sorted(k))))
     &                 .AND.((.not.cv_constr_mem).or.
     &                 (cv_proc_memused(cv_proc_sorted(k))+aux_mem.lt.
     &                 cv_proc_maxmem(cv_proc_sorted(k))))
     &                 .AND.    
     &                 (cv_layer_p2node(layernmb)%
     &                 t2_cand(j,cv_proc_sorted(k)).eq.0))
     &                 then
                     cv_proc_workload(cv_proc_sorted(k))=
     &                    cv_proc_workload(cv_proc_sorted(k))+aux_flop
                     cv_proc_memused(cv_proc_sorted(k))=
     &                    cv_proc_memused(cv_proc_sorted(k))+aux_mem
                     cv_layer_p2node(layernmb)%
     &                    t2_cand(j,cv_proc_sorted(k))
     &                    =inode
                     cv_layerworkload(cv_proc_sorted(k))=
     &                    cv_layerworkload(cv_proc_sorted(k))+aux_flop
                     cv_layermemused(cv_proc_sorted(k))=
     &                    cv_layermemused(cv_proc_sorted(k))+aux_mem
                     nmb_cand_needed=nmb_cand_needed-1
                     k=k+1
                  else
                     k=k+1
                     if(k.gt.cv_slavef) then
                        if(cv_lp.gt.0)
     &                       write(cv_lp,*)'Error in ',subname
                        goto 999
                     endif
                  end if
               end do
               if(nmb_cand_needed.gt.0) then
                  if(cv_lp.gt.0)
     &                 write(cv_lp,*)'Error in ',subname
                  goto 999
               endif
            end if
         end if
      end do
      do i=1,cv_layer_p2node(layernmb)%nmb_t2s
         nmb_cand_needed=
     &                  cv_layer_p2node(layernmb)%t2_cand(i,cv_slavef+1)
         candid= cv_layer_p2node(layernmb)%t2_cand(i,1:cv_slavef)
         cv_layer_p2node(layernmb)%t2_cand(i,1:cv_slavef)=-1
         k=0
         do j=1,cv_slavef
            if(candid(j).gt.0) then
               k=k+1
               cv_layer_p2node(layernmb)%t2_cand(i,k)=j-1
            end if
         end do
         if (k.ne.nmb_cand_needed) then
            if(cv_lp.gt.0)
     &           write(cv_lp,*)'Error in ',subname
            goto 999
         endif
      enddo
      do i=1,cv_slavef
         cv_layerworkload(i)=cv_layerworkload(i)-old_workload(i)
         cv_layermemused(i)=cv_layermemused(i)-old_memused(i)
      enddo
      istat=0
 999  continue
      DEALLOCATE(candid, sorted_nmb, sorted_costw, sorted_costm,
     & old_workload, old_memused)
      return
      end subroutine MUMPS_MAP_LAYER
      recursive subroutine MUMPS_MAPBELOW(inode,procnmb,
     &                                       procnode)
      integer,intent(in)::inode,procnmb
      integer,intent(inout)::procnode(:) 
      integer in
      procnode(inode)=procnmb
      if (cv_fils(inode).eq.0) return 
      in=cv_fils(inode)
      do while(in>0)
         procnode(in)=procnmb
         in=cv_fils(in)
      end do
      in=-in
      do while(in>0)
         call MUMPS_MAPBELOW(in,procnmb,procnode)
         in=cv_frere(in)
      end do
      return
      end subroutine MUMPS_MAPBELOW
      subroutine MUMPS_MAPSUBTREE(procnode)
      implicit none
      integer,intent(inout)::procnode(:) 
      integer i,inode,procnmb
      do i=cv_layerl0_start,cv_layerl0_end
         inode=cv_layerl0_array(i)
         if(inode.gt.0) then
            procnmb=procnode(inode)
            call MUMPS_MAPBELOW(inode,procnmb,procnode)
         endif
      enddo
      return
      end subroutine MUMPS_MAPSUBTREE
      subroutine MUMPS_POSTPROCESS_MEM()
      implicit none
      integer candid,inode,index,i,j,layernmb,master,nmbcand,swapper,
     &        totalnmb,node_of_master,node_of_candid,node_of_swapper
      DOUBLE PRECISION::mastermem,slavemem,maxmem
      logical swapthem,cand_better_master_arch,cand_better_swapper_arch
      intrinsic maxval,minval
      maxmem=maxval(cv_proc_memused(:))
      totalnmb=0
      do layernmb=cv_maxlayer,1,-1
         do i=1,cv_layer_p2node(layernmb)%nmb_t2s
            inode=cv_layer_p2node(layernmb)%t2_nodenumbers(i)
            master=cv_procnode(inode)
            if(ke69 .gt. 1) then
               allowed_nodes = .FALSE.
               call MUMPS_FIX_ACCEPTED_MASTER(layernmb,i)
               node_of_master = mem_distribmpi(master-1)
               if (node_of_master .lt. 0 ) then
                if(cv_mp.gt.0) write(cv_mp,*)'node_of_master_not found'
               endif
               node_of_swapper = node_of_master
            endif
            mastermem=cv_proc_memused(master)
            nmbcand=cv_layer_p2node(layernmb)%t2_cand(i,cv_slavef+1)
            swapper=master
            index=0
            do j=1,nmbcand
               candid=cv_layer_p2node(layernmb)%t2_cand(i,j)+1
               slavemem=cv_proc_memused(candid)
               if(ke69 .gt. 1) then
                  node_of_candid = mem_distribmpi(candid-1)
                  if (node_of_candid .lt. 0 ) then
                     if(cv_mp.gt.0) write(cv_mp,*)
     &               'node_of_candid_not found'
                  endif
               endif
               if(ke69 .le. 1) then
                  if((slavemem.lt.mastermem) .and.
     &                 (slavemem.lt.cv_proc_memused(swapper))) then
                     swapper=candid
                     index=j
                  endif
               else
                  cand_better_master_arch = (
     &                    (
     &                       (slavemem.lt.mastermem) .or.
     &                       (.not. allowed_nodes(node_of_master))
     &                    )
     &                    .and. allowed_nodes(node_of_candid)
     &                 )
                  cand_better_swapper_arch = (
     &                    (
     &                       (slavemem.lt.cv_proc_memused(swapper)) .or.
     &                       (.not. allowed_nodes(node_of_swapper))
     &                    )
     &                    .and. allowed_nodes(node_of_candid)
     &                 )
                  if(cand_better_master_arch .and.
     &                 cand_better_swapper_arch  ) then
                     swapper=candid
                     node_of_swapper = node_of_candid
                     index=j
                  endif
               endif
            enddo
            if(swapper.ne.master) then
               swapthem = .FALSE.
               if(0.75D0*mastermem.ge.cv_proc_memused(swapper))
     &            swapthem=.TRUE.
               if(mastermem.le.mastermem-cv_ncostm(inode)
     &                +cv_layer_p2node(layernmb)%t2_candcostm(i))
     &            swapthem=.FALSE.
               if(mastermem.le.cv_proc_memused(swapper)
     &                    +cv_ncostm(inode)
     &                    -cv_layer_p2node(layernmb)%t2_candcostm(i))
     &            swapthem=.FALSE.
               if(maxmem.le.mastermem-cv_ncostm(inode)
     &                +cv_layer_p2node(layernmb)%t2_candcostm(i))
     &            swapthem=.FALSE.
               if(maxmem.le.cv_proc_memused(swapper)+cv_ncostm(inode)
     &                -cv_layer_p2node(layernmb)%t2_candcostm(i))
     &            swapthem=.FALSE.
               if(ke69 .gt. 1) then
                  if (.not. allowed_nodes(node_of_master)) then
                     swapthem=.TRUE.
                  endif
               endif
               if(.NOT.swapthem) cycle
               cv_proc_workload(master)=cv_proc_workload(master)
     &                -cv_ncostw(inode)
     &                +cv_layer_p2node(layernmb)%t2_candcostw(i)
               cv_proc_memused(master)=cv_proc_memused(master)
     &                -cv_ncostm(inode)
     &                +cv_layer_p2node(layernmb)%t2_candcostm(i)
               cv_proc_workload(swapper)=cv_proc_workload(swapper)
     &                +cv_ncostw(inode)
     &                -cv_layer_p2node(layernmb)%t2_candcostw(i)
               cv_proc_memused(swapper)=cv_proc_memused(swapper)
     &                +cv_ncostm(inode)
     &                -cv_layer_p2node(layernmb)%t2_candcostm(i)
               cv_layer_p2node(layernmb)%t2_cand(i,index)=master-1
               cv_procnode(inode)=swapper
               maxmem=maxval(cv_proc_memused(:))
               totalnmb = totalnmb+1
            endif
         enddo
      enddo
      end subroutine MUMPS_POSTPROCESS_MEM
      subroutine MUMPS_PROCINIT(maxwork,maxmem,istat)
      implicit none
      DOUBLE PRECISION,intent(in),OPTIONAL::maxwork(cv_slavef),
     &                                      maxmem(cv_slavef)
      integer,intent(out)::istat
      integer i,allocok
      intrinsic huge
      DOUBLE PRECISION dummy
      character (len=48):: subname
      istat=-1
      subname='PROCINIT'
      if(present(maxwork)) then
         cv_constr_work=.TRUE.
      else
         cv_constr_work=.FALSE.
      end if
      if(present(maxmem)) then
         cv_constr_mem=.TRUE.
      else
         cv_constr_mem=.FALSE.
      end if
      allocate(cv_proc_workload(cv_slavef),
     &         cv_proc_maxwork(cv_slavef),
     &         cv_proc_memused(cv_slavef),
     &         cv_proc_maxmem(cv_slavef),
     &         cv_proc_sorted(cv_slavef),
     &         STAT=allocok)
      if (allocok.gt.0) then
         cv_info(1) = cv_error_memalloc
         cv_info(2) = 2*cv_slavef
         istat = cv_error_memalloc
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'memory allocation error in ',subname
         return
      end if
      allocate(work_per_proc(cv_slavef),id_son(cv_slavef),STAT=allocok)
      if (allocok.gt.0) then
         cv_info(1) = cv_error_memalloc
         cv_info(2) = 2*cv_slavef
         istat = cv_error_memalloc
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'memory allocation error in ',subname
         return
      end if
      do i=1,cv_slavef
         cv_proc_workload(i)=dble(0)
         if(cv_constr_work) then
            cv_proc_maxwork(i)=maxwork(i)
         else
            cv_proc_maxwork(i)=(huge(dummy))
         endif
         cv_proc_memused(i)=dble(0)
         if(cv_constr_mem) then
            cv_proc_maxmem(i)=maxmem(i)
         else
            cv_proc_maxmem(i)=(huge(dummy))
         endif
      end do
      do i=1, cv_slavef
        cv_proc_sorted(i)=i
      enddo
      istat=0
      return
      end subroutine MUMPS_PROCINIT
      recursive subroutine MUMPS_MOD_PROPMAP
     &                    (inode,ctr,istat)
      implicit none
      integer, intent(in)::inode,ctr
      integer, intent(inout)::istat
      integer::j,k,in,in1,ierr,son,nmb_procs_inode,nmb_sons_inode,
     &         current,i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: procs4son 
      INTEGER :: allocok
      character (len=48):: subname
      DOUBLE PRECISION :: relative_weight,costs_sons
      DOUBLE PRECISION :: loc_relax
      INTEGER :: depth
      logical force_cand
      DOUBLE PRECISION Y
      intrinsic random_number
      integer nmb_propmap_strict,share2,procsrest,current2
      integer k69onid
      INTEGER, ALLOCATABLE, DIMENSION(:) :: procs_inode
      LOGICAL UPDATE_CTR
      if (ctr.le.0) then
         istat = 0
         return
      endif
      istat= -1
      if(cv_frere(inode).eq.cv_n+1) return
      subname='MOD_PROPMAP'
      if(.NOT.associated(cv_prop_map(inode)%ind_proc)) return
      ALLOCATE(procs_inode(cv_slavef),
     &         procs4son(cv_size_ind_proc),stat=allocok)
      if (allocok.gt.0) then
         cv_info(1) = cv_error_memalloc
         cv_info(2) = cv_size_ind_proc + cv_slavef
         istat = cv_error_memalloc
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'memory allocation error in ',subname
         return
      end if
      procs_inode=-1
      nmb_procs_inode = 0
      do j=1,cv_slavef
         if( MUMPS_BIT_GET4PROC(inode,j))then
            nmb_procs_inode = nmb_procs_inode + 1
         endif
      end do
      i=0
      do j=1,cv_slavef
            if(ke69 .gt.1) then
               call MUMPS_GET_IDP1_PROC(j-1,
     &              k69onid,ierr)
            else
               k69onid = j
            endif
            if(MUMPS_BIT_GET4PROC(inode,k69onid))then
               i = i + 1
               procs_inode(i)=k69onid
            endif
      end do
      if(i.ne.nmb_procs_inode)then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error in ',subname
     &        ,subname
         goto 999
      endif
      if(nmb_procs_inode.eq.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error in ',subname
     &        ,subname
         goto 999
      end if
      if ((cv_nodelayer(inode).eq.0).AND.
     &         (cv_frere(inode).ne.cv_n+1)) then
         istat = 0
         goto 999
      endif
      nmb_sons_inode = 0
      costs_sons = dble(0)
      force_cand=(mod(cv_keep(24),2).eq.0)
      in = inode
      do while (cv_fils(in).gt.0)
         in=cv_fils(in)
      end do
      if (cv_fils(in).eq.0) then
         istat = 0
         goto 999
      endif
      in = -cv_fils(in)
      son=in 
      do while(in.gt.0)
         nmb_sons_inode = nmb_sons_inode + 1
         if(cv_tcostw(in).le.0.0D0) then
            if(cv_lp.gt.0)
     &           write(cv_lp,*)'Subtree costs for ',in,
     &                         ' should be positive in ',subname
            goto 999
         endif
         costs_sons = costs_sons + cv_tcostw(in)
         in=cv_frere(in)
      enddo
      if(costs_sons.le.0D0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error in ',subname
     &        ,subname
         goto 999
      endif
      depth= max(cv_mixed_strat_bound - ctr,0)
      if ((cv_keep(24).eq.16).OR.(cv_keep(24).eq.18)) then
         if(depth.ge.cv_mixed_strat_bound) then
            loc_relax = dble(1)
         else
            loc_relax =  dble(1) +
     &            max(dble(cv_keep(77))/dble(100), dble(0))
         endif
      else
         loc_relax = dble(1)
      endif
      in=son
      current = 1
      do while(in.gt.0)
         UPDATE_CTR = .TRUE.
         if( ( (nmb_sons_inode.ge.nmb_procs_inode).AND.
     &         (nmb_procs_inode.LT.4) )
     &       .OR. ( nmb_sons_inode.EQ.1 )
     &     ) then
            procs4son = cv_prop_map(inode)%ind_proc
            IF (nmb_sons_inode.EQ.1) UPDATE_CTR=.FALSE.
         else
            do k=1,cv_size_ind_proc
               do j=0,cv_bitsize_of_int-1
                  procs4son(k)=ibclr(procs4son(k),j)
               end do
            end do
            nmb_propmap_strict=0
            do k=1,cv_slavef
               if( MUMPS_BIT_GET4PROC(in,k)) then
                  nmb_propmap_strict=nmb_propmap_strict+1
                  call MUMPS_BIT_SET(procs4son,k,ierr)
               end if
            end do
            if(costs_sons.gt.0.0D0) then
               relative_weight=cv_tcostw(in)/costs_sons
            else
               relative_weight=0.0D0
            endif
            current = nmb_propmap_strict
            share2=
     &           max(0,nint(relative_weight*(loc_relax-dble(1))*
     &           dble(nmb_procs_inode)))
            procsrest=nmb_procs_inode - nmb_propmap_strict
            share2=min(share2,procsrest)
            CALL random_number(Y)
            current2=int(dble(Y)*dble(procsrest))
            k=1
            i=1
            do while((share2.gt.0).and.(i.le.2))
               do j=1,nmb_procs_inode
                  if(share2.le.0) exit
                  k69onid = procs_inode(j)
                  if(( MUMPS_BIT_GET4PROC(inode,k69onid)).AND.
     &                 (.NOT.MUMPS_BIT_GET(procs4son,k69onid))) then
                     if(k.ge.current2)then
                        call MUMPS_BIT_SET(procs4son,k69onid,ierr)
                        if(ierr.ne.0) then
                           if(cv_lp.gt.0)write(cv_lp,*)
     &                          'BIT_SET signalled error to',subname
                           istat = ierr
                           goto 999
                        end if
                        share2 = share2 - 1
                     endif
                     k=k+1
                  end if
               enddo
               i=i+1
            enddo
            if(share2.ne.0) then
               if(cv_lp.gt.0) write(cv_lp,*)
     &           'Error reported in ',subname
               goto 999
            end if
         end if
         ierr=0
         in1=in
         cv_prop_map(in1)%ind_proc=procs4son
         IF (UPDATE_CTR) THEN
          call MUMPS_MOD_PROPMAP(in1,ctr-1,ierr)
         ELSE
          call MUMPS_MOD_PROPMAP(in1,ctr,ierr)
         ENDIF
         if(ierr.ne.0) then
            if(cv_lp.gt.0) write(cv_lp,*)
     &           'Error reported in ',subname
            istat=ierr
            goto 999
         endif
         in=cv_frere(in)
      end do
      istat = 0
 999  continue
      DEALLOCATE(procs_inode,procs4son)
      return
      end subroutine MUMPS_MOD_PROPMAP
      recursive subroutine MUMPS_PROPMAP(inode,ctr,istat)
      implicit none
      integer, intent(in)::inode,ctr
      integer, intent(inout)::istat
      integer::j,k,in,in1,ierr,son,nmb_procs_inode,nmb_sons_inode,
     &           share,current,offset,
     &           in_tmp,nfront,npiv,ncb,
     &           keep48_loc,min_cand_needed
      integer, dimension(:), allocatable :: procs4son
      character (len=48):: subname
      DOUBLE PRECISION :: relative_weight,costs_sons, shtemp
      DOUBLE PRECISION :: costs_sons_real
      DOUBLE PRECISION :: PartofaProc
      LOGICAL          :: SkipSmallNodes
      PARAMETER (PartofaProc=0.01D0)
      DOUBLE PRECISION :: loc_relax
      INTEGER :: depth
      logical force_cand
      integer MUMPS_REG_GETKMAX, MUMPS_BLOC2_GET_NSLAVESMIN
      external MUMPS_REG_GETKMAX, MUMPS_BLOC2_GET_NSLAVESMIN
      DOUBLE PRECISION Y
      intrinsic random_number
      integer nmb_propmap_strict,share2,procsrest,current2
      integer k69onid,nb_free_procs,local_son_indice,nb_procs_for_sons,
     &     ptr_upper_ro_procs
      INTEGER :: allocok
      logical upper_round_off,are_sons_treated
      DOUBLE PRECISION tmp_cost
      if (ctr.le.0) then
         istat = 0
         return
      endif
      istat= -1
      if(cv_frere(inode).eq.cv_n+1) return
      subname='PROPMAP'
      nmb_procs_inode = 0
      do j=1,cv_slavef
         if( MUMPS_BIT_GET4PROC(inode,j))
     &        nmb_procs_inode = nmb_procs_inode + 1
      end do
      if(nmb_procs_inode.eq.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error in ',subname
     &        ,subname
         return
      end if
      if ((cv_nodelayer(inode).eq.0).AND.
     &         (cv_frere(inode).ne.cv_n+1)) then
         istat = 0
         return
      endif
      ptr_upper_ro_procs=1
      work_per_proc(1:cv_slavef)=0.0D0
      id_son(1:cv_slavef)=0
      nmb_sons_inode = 0
      costs_sons = dble(0)
      force_cand=(mod(cv_keep(24),2).eq.0)
      min_cand_needed=0
      in = inode
      do while (cv_fils(in).gt.0)
         in=cv_fils(in)
      end do
      if (cv_fils(in).eq.0) then
         istat = 0
         return
      endif
      in = -cv_fils(in)
      son=in 
      do while(in.gt.0)
         nmb_sons_inode = nmb_sons_inode + 1
         if(cv_tcostw(in).le.0.0D0) then
            if(cv_lp.gt.0)
     &           write(cv_lp,*)'Subtree costs for ',in,
     &                         ' should be positive in ',subname
            return
         endif
         costs_sons = costs_sons + cv_tcostw(in)
         in=cv_frere(in)
      enddo
      costs_sons_real = costs_sons
      SkipSmallNodes = .true.
      IF (costs_sons_real.gt.0.0D0) then
       in = son
       do while (in.gt.0) 
          relative_weight=cv_tcostw(in)/costs_sons_real
          shtemp = relative_weight*dble(nmb_procs_inode)
          IF (shtemp.lt.PartofaProc) THEN
            costs_sons = costs_sons - cv_tcostw(in)
          ENDIF
          in=cv_frere(in)
       enddo
       IF (costs_sons.LT. PartofaProc*costs_sons_real) THEN
         costs_sons = costs_sons_real
         SkipSmallNodes = .false.
       ENDIF
      ENDIF
      if(costs_sons.le.0.0D0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error in ',subname
     &        ,subname
         return
      endif
      if(cv_relax.le.0.0D0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error in ',subname,'. Wrong cv_relax'
         return
      endif
      ALLOCATE(procs4son(cv_size_ind_proc),stat=allocok)
      IF (allocok .GT. 0) THEN
         cv_info(1) = cv_error_memalloc
         cv_info(2) = cv_size_ind_proc
         istat = cv_error_memalloc
         if(cv_lp.gt.0)
     &           write(cv_lp,*)
     &           'Memory allocation error in ',subname
         return
      ENDIF
      depth= max(cv_n - ctr,0)
      if(cv_keep(24).eq.8) then
         loc_relax = cv_relax
      elseif ((cv_keep(24).eq.16).OR.(cv_keep(24).eq.18)) then
         loc_relax = cv_relax
      elseif (cv_keep(24).eq.10) then
         loc_relax = cv_relax
      elseif ((cv_keep(24).eq.12).OR.(cv_keep(24).eq.14)) then
         if(depth.ge.cv_mixed_strat_bound) then
            loc_relax = cv_relax
         else
            loc_relax =  cv_relax +
     &          max(dble(cv_keep(77))/dble(100), dble(0))
         endif
      endif
      in=son
      current = 1
      local_son_indice=1
      nb_procs_for_sons=0
      upper_round_off=.FALSE.
      are_sons_treated=.TRUE.
      do while(in.gt.0)
         if( (nmb_sons_inode.ge.nmb_procs_inode).AND.
     &       (nmb_procs_inode.LT.4) ) then
            procs4son = cv_prop_map(inode)%ind_proc
            are_sons_treated=.FALSE.
            nb_procs_for_sons=nmb_procs_inode
            nmb_propmap_strict=nmb_procs_inode
         elseif(nmb_procs_inode .LE. cv_keep(83)) then
            procs4son = cv_prop_map(inode)%ind_proc
            are_sons_treated=.FALSE.
            nb_procs_for_sons=nmb_procs_inode
            nmb_propmap_strict=nmb_procs_inode
         else
            do k=1,cv_size_ind_proc
               do j=0,cv_bitsize_of_int-1
                  procs4son(k)=ibclr(procs4son(k),j)
               end do
            end do
            if(costs_sons.gt.0.0D0) then
               relative_weight=cv_tcostw(in)/costs_sons
            else
               relative_weight=dble(0)
            endif
            shtemp = relative_weight*dble(nmb_procs_inode)
            IF ( (shtemp.LT.PartofaProc)
     &             .AND. ( SkipSmallNodes ) )  THEN
             share = 1  
             do j=current,cv_slavef
              if(ke69 .gt.1) then
               call MUMPS_GET_IDP1_PROC(j-1,k69onid,ierr)
              else
               k69onid = j
              endif 
              if( MUMPS_BIT_GET4PROC(inode,k69onid)) then
                  call MUMPS_BIT_SET(procs4son,k69onid,ierr)
                  if(ierr.ne.0) then
                     if(cv_lp.gt.0)write(cv_lp,*)
     &               'BIT_SET signalled error to',subname
                     istat = ierr
                     goto 999
                  end if
                  share = share -1
                  exit 
              endif
             enddo
             if (share.gt.0) then
               do j=1,current-1
                if(ke69 .gt.1) then
                 call MUMPS_GET_IDP1_PROC(j-1,k69onid,ierr)
                else
                 k69onid = j
                endif 
                if( MUMPS_BIT_GET4PROC(inode,k69onid)) then
                    call MUMPS_BIT_SET(procs4son,k69onid,ierr)
                    if(ierr.ne.0) then
                       if(cv_lp.gt.0)write(cv_lp,*)
     &               'BIT_SET signalled error to',subname
                       istat = ierr
                       goto 999
                    end if
                    share = share -1
                    exit 
                endif
               enddo
              endif
              if(share.ne.0) then
                 if(cv_lp.gt.0) write(cv_lp,*)
     &              'Error reported in ',subname
                 goto 999
              end if
              if(.NOT.associated(cv_prop_map(in)%ind_proc)) then
                 call MUMPS_PROPMAP_INIT(in,ierr)
              if(ierr.ne.0) then
               if(cv_lp.gt.0)
     &              write(cv_lp,*)'PROPMAP_INIT signalled error to'
     &              ,subname
               istat = ierr
               goto 999
              end if
             endif
             current  = j 
             cv_prop_map(in)%ind_proc = procs4son
             in = cv_frere(in)
             cycle   
            ENDIF
            share  = max(1,nint(shtemp))
            if (dble(share).ge.shtemp) then
               upper_round_off=.TRUE.
            else
               upper_round_off = .FALSE.
            endif
            share=min(share,nmb_procs_inode)
            nmb_propmap_strict=share
            nb_procs_for_sons=nb_procs_for_sons+nmb_propmap_strict
            offset=1
            do j=current,cv_slavef
               if(ke69 .gt.1) then
                  call MUMPS_GET_IDP1_PROC(j-1,k69onid,ierr)
               else
                  k69onid = j
               endif
               if( MUMPS_BIT_GET4PROC(inode,k69onid)) then
                  call MUMPS_BIT_SET(procs4son,k69onid,ierr)
                  if(ierr.ne.0) then
                     if(cv_lp.gt.0)write(cv_lp,*)
     &               'BIT_SET signalled error to',subname
                     istat = ierr
                     goto 999
                  end if
                  share = share-1
                  if(share.le.0) then
                     current = j + offset
                     if(current.gt.cv_slavef) current = 1
                     exit
                  end if
               end if
            end do
            if(share.gt.0) then
               do j=1,current-1
                  if(ke69 .gt.1) then
                     call MUMPS_GET_IDP1_PROC(j-1,k69onid,ierr)
                  else
                     k69onid = j
                  endif
                  if( MUMPS_BIT_GET4PROC(inode,k69onid)) then
                     call MUMPS_BIT_SET(procs4son,k69onid,ierr)
                     if(ierr.ne.0) then
                        if(cv_lp.gt.0)write(cv_lp,*)
     &                  'BIT_SET signalled error to',subname
                        istat = ierr
                        goto 999
                     end if
                     share = share-1
                     if(share.le.0) then
                        current = j + offset
                        if(current.gt.cv_slavef) current = 1
                        exit
                     end if
                  end if
               end do
            endif
            if(share.ne.0) then
               if(cv_lp.gt.0) write(cv_lp,*)
     &              'Error reported in ',subname
               goto 999
            end if
            if(.not.upper_round_off)then
               if(local_son_indice.lt.cv_slavef)then
                  id_son(local_son_indice)=in
                  work_per_proc(local_son_indice)=cv_tcostw(in)/
     &                 dble(nmb_propmap_strict)
                  local_son_indice=local_son_indice+1
                  if(local_son_indice.eq.cv_slavef)then
                     CALL MUMPS_SORT_MSORT(ierr,cv_slavef,id_son,
     &                    work_per_proc)
                     if(ierr.ne.0) then
                       if(cv_lp.gt.0)
     &                   write(cv_lp,*)
     &                  'Error reported by MUMPS_SORT_MSORT in ',subname
                         istat = ierr
                         goto 999
                       endif
                  endif
               else
                  current2=cv_slavef
                  tmp_cost=cv_tcostw(in)/dble(nmb_propmap_strict)
                  do while(current2.ge.1)
                     if(tmp_cost.lt.work_per_proc(current2))exit
                     current2=current2-1
                  enddo
                  if(current2.ne.cv_slavef)then
                     if(current2.eq.0)then
                        current2=1
                     endif
                     do j=cv_slavef-1,current2,-1
                        id_son(j+1)=id_son(j)
                        work_per_proc(j+1)=work_per_proc(j)
                     enddo
                     id_son(current2)=in
                     work_per_proc(current2)=tmp_cost
                  endif
               endif
            endif
            upper_round_off=.FALSE.
         endif
         if(.NOT.associated(cv_prop_map(in)%ind_proc)) then
            call MUMPS_PROPMAP_INIT(in,ierr)
            if(ierr.ne.0) then
               if(cv_lp.gt.0)
     &              write(cv_lp,*)'PROPMAP_INIT signalled error to'
     &              ,subname
               istat = ierr
               goto 999
            end if
         endif
         cv_prop_map(in)%ind_proc = procs4son
         in=cv_frere(in)
      end do
      if(are_sons_treated)then
         if(nb_procs_for_sons.ne.nmb_procs_inode)then
            do j=1,nmb_procs_inode-nb_procs_for_sons
               procs4son=cv_prop_map(id_son(j))%ind_proc
               do while(current.le.cv_slavef)
                  if(ke69 .gt.1) then
                     call MUMPS_GET_IDP1_PROC(current-1,k69onid,ierr)
                  else
                     k69onid = current
                  endif
                  if(.NOT.MUMPS_BIT_GET4PROC(inode,k69onid)) then
                     current=current+1
                  else
                     exit  
                  endif
               enddo
               call MUMPS_BIT_SET(procs4son,k69onid,ierr)
               cv_prop_map(id_son(j))%ind_proc=procs4son
            enddo
            ptr_upper_ro_procs=min(j,nmb_procs_inode-nb_procs_for_sons)
         endif
      endif
      in=son
      current = 1
      do while(in.gt.0)
         if( (nmb_sons_inode.ge.nmb_procs_inode).AND.
     &       (nmb_procs_inode.LT.4) ) then
            procs4son = cv_prop_map(inode)%ind_proc
         elseif(nmb_procs_inode .LE. cv_keep(83)) then
            procs4son = cv_prop_map(inode)%ind_proc
         else
            procs4son = cv_prop_map(in)%ind_proc
            in_tmp=in
            nfront=cv_nfsiz(in_tmp)
            npiv=0
            in_tmp=in_tmp
            do while(in_tmp.gt.0)
              if (cv_BLKON) then
                npiv = npiv + cv_SIZEOFBLOCKS(in_tmp)
              else
               npiv=npiv+1
              endif
              in_tmp=cv_fils(in_tmp)
            end do
            ncb=nfront-npiv
            if (force_cand) then
               if (cv_keep(50) ==  0) then
                  keep48_loc=0
               else
                  keep48_loc=3
               endif
               if (cv_keep(48).EQ.5) keep48_loc = 5
               min_cand_needed=
     &              MUMPS_BLOC2_GET_NSLAVESMIN
     &              (cv_slavef, keep48_loc,cv_keep8(21),
     &              cv_keep(50),
     &              nfront,ncb,
     &              cv_keep(375), cv_keep(119))
               min_cand_needed=min(cv_slavef,min_cand_needed+1)
            else
               min_cand_needed = 0
            endif
            min_cand_needed = max(min_cand_needed, cv_keep(91))
            if(costs_sons.gt.0.0D0) then
               relative_weight=cv_tcostw(in)/costs_sons
            else
               relative_weight=dble(0)
            endif
            nmb_propmap_strict=0
            do k=1,cv_slavef
               if( MUMPS_BIT_GET(procs4son,k)) then
                  nmb_propmap_strict=nmb_propmap_strict+1
               end if
            end do
            offset=1
            share2=
     &          max(0,nint(relative_weight*(loc_relax-dble(1))*
     &                                   dble(nmb_procs_inode)))
            share2 = max(share2, min_cand_needed -nmb_propmap_strict,
     &                   (cv_keep(83)/2) - nmb_propmap_strict)
            procsrest=nmb_procs_inode - nmb_propmap_strict
            share2=min(share2,procsrest)
            share2 = 0 
            CALL random_number(Y)
            current2     =int(dble(Y)*dble(procsrest))
            nb_free_procs=1
            do j=1,cv_slavef
               if(share2.le.0) exit
               if(ke69 .gt.1) then
                     call MUMPS_GET_IDP1_PROC(j-1,k69onid,ierr)
                  else
                     k69onid = j
                  endif
               if(( MUMPS_BIT_GET4PROC(inode,k69onid)).AND.
     &           (.NOT.MUMPS_BIT_GET(procs4son,k69onid))) then
                  if(nb_free_procs.ge.current2)then
                     call MUMPS_BIT_SET(procs4son,k69onid,ierr)
                     if(ierr.ne.0) then
                        if(cv_lp.gt.0)write(cv_lp,*)
     &                       'BIT_SET signalled error to',subname
                        istat = ierr
                        goto 999
                     end if
                     share2 = share2 - 1
                  endif
                  nb_free_procs=nb_free_procs+1
               end if
            end do
            if(share2.gt.0) then
               do j=1,cv_slavef
                  if(share2.le.0) exit
                  if(ke69 .gt.1) then
                     call MUMPS_GET_IDP1_PROC(j-1,k69onid,ierr)
                  else
                     k69onid = j
                  endif
                  if(( MUMPS_BIT_GET4PROC(inode,k69onid)).AND.
     &              (.NOT.MUMPS_BIT_GET(procs4son,k69onid))) then
                        call MUMPS_BIT_SET(procs4son,k69onid,ierr)
                        if(ierr.ne.0) then
                           if(cv_lp.gt.0)write(cv_lp,*)
     &                          'BIT_SET signalled error to',subname
                           istat = ierr
                           goto 999
                        end if
                        share2 = share2 - 1
                  end if
               end do
            endif
            if(share2.ne.0) then
               if(cv_lp.gt.0) write(cv_lp,*)
     &              'Error reported in ',subname
               goto 999
            end if
         endif
         ierr=0
         in1=in
         cv_prop_map(in1)%ind_proc = procs4son
         call MUMPS_PROPMAP(in1,ctr-1,ierr)
         if(ierr.ne.0) then
            if(cv_lp.gt.0) write(cv_lp,*)
     &         'Error reported in ',subname
            istat=ierr
            goto 999
         endif
         in=cv_frere(in)
      end do
      istat = 0
 999  CONTINUE
      DEALLOCATE(procs4son)
      return
      end subroutine MUMPS_PROPMAP
      subroutine MUMPS_PROPMAP_INIT(inode,istat)
      implicit none
      integer, intent(in)::inode
      integer, intent(out)::istat
      integer j,k,allocok
      character (len=48):: subname
      istat = -1
      if(cv_frere(inode).eq.cv_n+1) return
      subname='PROPMAP_INIT'
      if(.not.associated(
     &     cv_prop_map(inode)%ind_proc)) then
         allocate(cv_prop_map(inode)%ind_proc
     &        (cv_size_ind_proc),STAT=allocok)
         if (allocok.gt.0) then
            cv_info(1) = cv_error_memalloc
            cv_info(2) = cv_size_ind_proc
            istat = cv_error_memalloc
            if(cv_lp.gt.0)
     &           write(cv_lp,*)
     &           'memory allocation error in ',subname
            return
         end if
      end if
      do k=1,cv_size_ind_proc
         do j=0,cv_bitsize_of_int-1
            cv_prop_map(inode)%ind_proc(k)=
     &         ibclr(cv_prop_map(inode)%ind_proc(k),j)
         end do
      end do
      istat = 0
      return
      end subroutine MUMPS_PROPMAP_INIT
      subroutine MUMPS_PROPMAP_TERM(inode,istat)
      integer,intent(in)::inode
      integer,intent(out)::istat
      integer ierr
      character (len=48):: subname
      subname='PROPMAP_TERM'
      istat =-1
      if(associated(cv_prop_map(inode)%ind_proc)) then
        deallocate(cv_prop_map(inode)%ind_proc, STAT=ierr)
        if(ierr.ne.0) then
           if(cv_lp.gt.0)
     &         write(cv_lp,*)'Memory deallocation error in ', subname
               istat = cv_error_memdeloc
               return
        endif
        nullify(cv_prop_map(inode)%ind_proc)
      end if
      istat =0
      return
      end subroutine MUMPS_PROPMAP_TERM
      subroutine MUMPS_PROPMAP4SPLIT(inode,ifather,istat)
      implicit none
      integer,intent(in)::inode,ifather
      integer,intent(out)::istat
      character (len=48):: subname
      istat= -1
      subname='PROPMAP4SPLIT'
      if((cv_frere(inode).eq.cv_n+1).OR.(cv_frere(ifather).eq.cv_n+1)
     &     .OR.(.NOT.associated(cv_prop_map(inode)%ind_proc))) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'tototo signalled error to'
     &        ,subname
         return
      endif
      if(.NOT.associated(cv_prop_map(ifather)%ind_proc)) then
         call MUMPS_PROPMAP_INIT(ifather,ierr)
         if(ierr.ne.0) then
            if(cv_lp.gt.0)
     &           write(cv_lp,*)'PROPMAP_INIT signalled error to '
     &           ,subname
            istat = ierr
            return
         end if
      endif
      cv_prop_map(ifather)%ind_proc =
     &                    cv_prop_map(inode)%ind_proc
      istat=0
      return
      end subroutine MUMPS_PROPMAP4SPLIT
      subroutine MUMPS_ROOTLIST(istat)
      implicit none
      integer,intent(out)::istat
      integer i,allocok
      character (len=48):: subname
      istat=-1
      subname='ROOTLIST'
      allocate(cv_layerl0_array(cv_maxnsteps),
     &         cv_layerl0_sorted_costw(cv_maxnsteps),STAT=allocok)
      if (allocok.gt.0) then
         cv_info(1) = cv_error_memalloc
         cv_info(2) = 12*cv_maxnsteps
         istat = cv_error_memalloc
         if(cv_lp.gt.0)
     &        write(cv_lp,*)
     &        'memory allocation error in ',subname
         return
      end if
      do i=1,cv_maxnsteps
         cv_layerl0_sorted_costw(i)=dble(0)
         cv_layerl0_array(i)=0
      end do
      cv_layerl0_start        = 0
      cv_layerl0_end          = 0
      layerL0_endforarrangeL0 = 0
      if ((.NOT.associated(cv_tcostw)).OR.(.NOT.associated(cv_tcostm)))
     &   then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error:tcost must be allocated in ',subname
         return
      end if
      cv_nbsa=0
      do i=1,cv_n
         if (cv_frere(i).eq.0) then
            cv_layerl0_start=1
            cv_layerl0_end=cv_layerl0_end+1
            IF (cv_tcostw(i).GT.mincostw)
     &           layerL0_endforarrangeL0 = layerL0_endforarrangeL0+1
            cv_layerl0_array(cv_layerl0_end)=i
            cv_layerl0_sorted_costw(cv_layerl0_end)=cv_tcostw(i)
            cv_costw_layer0=cv_costw_layer0 + cv_tcostw(i)
            cv_costm_layer0=cv_costm_layer0 + cv_tcostm(i)
            cv_nbsa=cv_nbsa+1
         end if
      end do
      if(cv_nbsa.eq.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Error:no root nodes in ',subname
         return
      end if
      call MUMPS_SORT_MSORT(ierr,cv_layerl0_end-cv_layerl0_start+1,
     &            cv_layerl0_array(cv_layerl0_start:cv_layerl0_end),
     &     cv_layerl0_sorted_costw(cv_layerl0_start:cv_layerl0_end))
      IF (ierr .ne.0) then
        if(cv_lp.gt.0)
     &    write(cv_lp,*)
     &    'Error reported by MUMPS_SORT_MSORT in ',subname
        istat = ierr
        return
      ENDIF
      cv_costw_total=cv_costw_layer0
      cv_costm_total=cv_costm_layer0
      istat=0
      return
      end subroutine MUMPS_ROOTLIST
      subroutine MUMPS_SELECT_TYPE3(istat)
      implicit none
      integer,intent(out)::istat
      character (len=48):: subname
      subname='SELECT_TYPE3'
      CALL MUMPS_SELECT_K38K20(cv_n, slavef, cv_mp, cv_icntl(13),
     &     cv_keep(1), cv_frere(1), cv_nfsiz(1), istat)
      IF (istat .NE. 0) THEN
            if(cv_lp.gt.0)
     &           write(cv_lp,*)
     &           'Error: Can''t select type 3 node in ',subname
      ELSE IF (cv_keep(38) .ne. 0) then
        IF(cv_nodelayer(cv_keep(38)).eq.0.and.
     &                            (cv_keep(60).EQ.0)) then
          cv_keep(38)=0
        ELSE
          cv_nodetype(cv_keep(38))=3
        ENDIF
      ENDIF
      RETURN
      end subroutine MUMPS_SELECT_TYPE3
      subroutine MUMPS_SETUP_CAND(istat)
      integer,intent(out):: istat
      integer            :: i,dummy,layernmb,allocok
      integer            :: montype, nbcand, inode
      character (len=48) :: subname
      istat=-1
      subname='SETUP_CAND'
      cv_nb_niv2=0
      do i=1,cv_n
         if(MUMPS_IS_NODE_OF_TYPE2(i)) cv_nb_niv2=cv_nb_niv2+1
      end do
      cv_keep(56)=cv_nb_niv2
      nullify(cv_par2_nodes,cv_cand)
      if(cv_nb_niv2.GT.0) then
         allocate(cv_par2_nodes(cv_nb_niv2),
     &        cv_cand(cv_nb_niv2,cv_slavef+1),STAT=allocok)
         if (allocok.gt.0) then
            cv_info(1) = cv_error_memalloc
            cv_info(2) = cv_nb_niv2*(cv_slavef+2)
            istat = cv_error_memalloc
            if(cv_lp.gt.0)
     &           write(cv_lp,*)
     &           'memory allocation error in ',subname
            return
         end if
         cv_par2_nodes=0
         cv_cand(:,:)=0
         dummy=1
         do layernmb=1,cv_maxlayer
            do i=1,cv_layer_p2node(layernmb)%nmb_t2s
               inode = cv_layer_p2node(layernmb)%t2_nodenumbers(i)
               cv_par2_nodes(dummy)= inode
               nbcand = cv_layer_p2node(layernmb)%t2_cand(i,cv_slavef+1)
               cv_cand(dummy,:)=cv_layer_p2node(layernmb)%t2_cand(i,:)
               montype= cv_nodetype(inode)
               if (montype.eq.tsplit_beg) then 
                  CALL MUMPS_SETUP_CAND_CHAIN(cv_n, cv_nb_niv2, 
     &                 cv_frere(1), cv_nodetype(1),
     &                 cv_par2_nodes(1), cv_procnode(1), cv_cand(1,1),
     &                 inode,
     &                 slavef, dummy, nbcand, istat)
               endif
               dummy=dummy+1
            enddo
         enddo
         if(dummy.ne.cv_nb_niv2+1) then
            if(cv_lp.gt.0)
     &           write(cv_lp,*)'Error in ',subname,
     &           ' : dummy =',dummy,'nbniv2 =',cv_nb_niv2
            return
         endif
      endif
      istat=0
      return
      end subroutine MUMPS_SETUP_CAND
      subroutine MUMPS_SORTPROCS(map_strat,workload,memused,
     &                       inode,istat)
      implicit none
      integer,intent(in)::map_strat
      DOUBLE PRECISION,dimension(:),intent(in)::workload, memused
      integer, optional::inode,istat
      integer i,j,aux_int,nmb_procs,pos
      character (len=48):: subname
      logical enforce_prefsort
      logical use_propmap
      logical,SAVE::init1 = .FALSE.
      logical,SAVE::init2 = .FALSE.
      subname='SORTPROCS'
      enforce_prefsort=.TRUE.
      use_propmap=present(inode)
      if(present(istat))istat=-1
      if((map_strat.ne.cv_equilib_flops).and.
     &   (map_strat.ne.cv_equilib_mem)) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'error in ',subname
         return
      endif
      i=0 
      do i = 1, cv_slavef
        cv_proc_sorted(i)=i
      enddo
      if (.not.present(inode)) then
         if(.NOT.init1) then
            init1=.TRUE.
         end if
         do i=1,cv_slavef-1
            do j=i+1,cv_slavef
               if(((workload(cv_proc_sorted(j)).lt.
     &                     workload(cv_proc_sorted(i))).AND.
     &             (map_strat.eq.cv_equilib_flops))
     &            .OR.
     &            ((memused(cv_proc_sorted(j)).lt.
     &                     memused(cv_proc_sorted(i))).AND.
     &             (map_strat.eq.cv_equilib_mem)))then
                  aux_int=cv_proc_sorted(j)
                  cv_proc_sorted(j)=cv_proc_sorted(i)
                  cv_proc_sorted(i)=aux_int
               end if
            end do
         end do
      else if(present(inode)) then
         if (use_propmap) then
            if(.NOT.init2) then
               init2=.TRUE.
            end if
            nmb_procs=0
            do pos=1,cv_slavef
               if( MUMPS_BIT_GET4PROC(inode,pos)) then
                  if (pos.le.nmb_procs) then
                     exit
                  else
                     nmb_procs=nmb_procs+1
                     aux_int=cv_proc_sorted(pos)
                     cv_proc_sorted(pos)=
     &                    cv_proc_sorted(nmb_procs)
                     cv_proc_sorted(nmb_procs)=aux_int
                     cycle
                  end if
               end if
            end do
         end if
         do i=1,nmb_procs-1
            do j=i+1,nmb_procs
               if(((workload(cv_proc_sorted(j)).lt.
     &                     workload(cv_proc_sorted(i))).AND.
     &             (map_strat.eq.cv_equilib_flops))
     &            .OR.
     &            ((memused(cv_proc_sorted(j)).lt.
     &                     memused(cv_proc_sorted(i))).AND.
     &             (map_strat.eq.cv_equilib_mem)))then
                  aux_int=cv_proc_sorted(j)
                  cv_proc_sorted(j)=cv_proc_sorted(i)
                  cv_proc_sorted(i)=aux_int
               end if
            end do
         end do
         do i=nmb_procs+1,cv_slavef-1
            do j=i+1,cv_slavef
               if(((workload(cv_proc_sorted(j)).lt.
     &                     workload(cv_proc_sorted(i))).AND.
     &             (map_strat.eq.cv_equilib_flops))
     &            .OR.
     &            ((memused(cv_proc_sorted(j)).lt.
     &                     memused(cv_proc_sorted(i))).AND.
     &             (map_strat.eq.cv_equilib_mem)))then
                  aux_int=cv_proc_sorted(j)
                  cv_proc_sorted(j)=cv_proc_sorted(i)
                  cv_proc_sorted(i)=aux_int
               end if
            end do
         end do
         if(.NOT.enforce_prefsort) then
            if(((2.0D0*workload(cv_proc_sorted(nmb_procs+1)).lt.
     &           workload(cv_proc_sorted(1))).AND.
     &          (map_strat.eq.cv_equilib_flops))
     &         .OR.
     &         ((2.0D0*memused(cv_proc_sorted(nmb_procs+1)).lt.
     &           memused(cv_proc_sorted(1))).AND.
     &          (map_strat.eq.cv_equilib_mem)))then
               do i=1,cv_slavef-1
                  do j=i+1,cv_slavef
                     if(((workload(cv_proc_sorted(j)).lt.
     &                    workload(cv_proc_sorted(i))).AND.
     &                    (map_strat.eq.cv_equilib_flops))
     &                    .OR.
     &                    ((memused(cv_proc_sorted(j)).lt.
     &                    memused(cv_proc_sorted(i))).AND.
     &                    (map_strat.eq.cv_equilib_mem)))then
                        aux_int=cv_proc_sorted(j)
                        cv_proc_sorted(j)=cv_proc_sorted(i)
                        cv_proc_sorted(i)=aux_int
                     end if
                  end do
               end do
            endif
         end if
      endif
      if(present(istat))istat=0
      return
      end subroutine MUMPS_SORTPROCS
      subroutine MUMPS_STORE_GLOBALS(ne,nfsiz,frere,fils,keep,KEEP8,
     &                                info,procnode,ssarbr,nbsa)
      implicit none
      integer,dimension(cv_n),intent(inout)::ne,nfsiz,frere,fils,
     &                                     procnode,ssarbr
      integer, intent(inout):: keep(500),info(80),nbsa
      INTEGER(8) KEEP8(150)
      ne=cv_ne
      nfsiz=cv_nfsiz
      frere=cv_frere
      fils=cv_fils
      keep(2) =cv_keep(2) 
      keep(20)=cv_keep(20)
      keep(28)=cv_nsteps
      keep(38)=cv_keep(38)
      keep(56)=cv_keep(56)
      keep(61)=cv_keep(61)
      info(5)=cv_info(5)  
      info(6)=cv_nsteps
      procnode=cv_procnode
      ssarbr=cv_ssarbr
      nbsa=cv_nbsa
      end subroutine MUMPS_STORE_GLOBALS
      subroutine MUMPS_TERMGLOB(istat)
      implicit none
      integer,intent(out)::istat
      integer i,ierr,layernmb
      character (len=48):: subname
      istat=-1
      subname='TERMGLOB'
      nullify(cv_frere,cv_fils,cv_nfsiz,cv_ne,cv_keep,cv_keep8,
     &        cv_icntl,cv_info,cv_procnode,cv_ssarbr)
      deallocate(cv_proc_workload,cv_proc_maxwork,cv_proc_memused,
     &    cv_proc_maxmem,cv_nodetype,
     &    cv_nodelayer,cv_proc_sorted,
     &    cv_ncostw,cv_ncostm,
     &    cv_layerworkload,cv_layermemused,
     &    STAT=ierr)
      if(ierr.ne.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Memory deallocation error in ',subname
         istat = cv_error_memdeloc
         return
      end if
      deallocate(work_per_proc,id_son,STAT=ierr)
      if(ierr.ne.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Memory deallocation error in ',subname
         istat = cv_error_memdeloc
         return
      end if
      do layernmb=1,cv_maxlayer
         if(cv_layer_p2node(layernmb)%nmb_t2s.gt.0) then
            deallocate(cv_layer_p2node(layernmb)%t2_nodenumbers,
     &                 cv_layer_p2node(layernmb)%t2_cand,
     &                 cv_layer_p2node(layernmb)%t2_candcostw,
     &                 cv_layer_p2node(layernmb)%t2_candcostm,
     &                 STAT=ierr)
            if(ierr.ne.0) then
               if(cv_lp.gt.0)
     &         write(cv_lp,*)'Memory deallocation error in ',
     &                        subname
               istat = cv_error_memdeloc
               return
            end if
         endif
      enddo
      if(associated(cv_layer_p2node)) then
         deallocate(cv_layer_p2node,STAT=ierr)
         if(ierr.ne.0) then
            if(cv_lp.gt.0)
     &      write(cv_lp,*)'Memory deallocation error in ',subname
            istat = cv_error_memdeloc
            return
         end if
      end if
      do i=1,cv_n
         call MUMPS_PROPMAP_TERM(i,ierr)
         if(ierr.ne.0) then
            if(cv_lp.gt.0)
     &           write(cv_lp,*)'PROPMAP_TERM signalled error in ',
     &           subname
            istat = ierr
            return
         end if
      end do
      if(associated(cv_prop_map))deallocate(cv_prop_map,STAT=ierr)
      if(ierr.ne.0) then
         if(cv_lp.gt.0)
     &        write(cv_lp,*)'Memory deallocation error in ',subname
         istat = cv_error_memdeloc
         return
      end if
      istat=0
      return
      end subroutine MUMPS_TERMGLOB
      recursive subroutine MUMPS_TREECOSTS(pos)
      implicit none
      integer,intent(in)::pos
      integer i,nfront,npiv,nextpos
      if ((.NOT.associated(cv_tcostw)).OR.(.NOT.associated(cv_tcostm)))
     &     then
         call MUMPS_ABORT()
      end if
      nfront=cv_nfsiz(pos)
      npiv=1
      nextpos=cv_fils(pos)
      do while (nextpos.gt.0)
         if (cv_BLKON) then
           npiv = npiv + cv_SIZEOFBLOCKS(nextpos)
         else
           npiv=npiv+1
         endif
         nextpos=cv_fils(nextpos)
      end do
      call MUMPS_CALCNODECOSTS(npiv,nfront,
     &     cv_ncostw(pos), cv_ncostm(pos))
      cv_tcostw(pos)=cv_ncostw(pos)
      cv_tcostm(pos)=cv_ncostm(pos)
      if (cv_ne(pos).ne.0) then
         nextpos=cv_fils(pos)
         do while(nextpos.gt.0)
            nextpos=cv_fils(nextpos)
         end do
         nextpos=-nextpos
         do i=1,cv_ne(pos)
            cv_depth(nextpos)=cv_depth(pos)+1
            call MUMPS_TREECOSTS(nextpos)
            cv_tcostw(pos)=cv_tcostw(pos)+cv_tcostw(nextpos)
            cv_tcostm(pos)=cv_tcostm(pos)+cv_tcostm(nextpos)
            nextpos=cv_frere(nextpos)
         end do
      endif
      return
      end subroutine MUMPS_TREECOSTS
      recursive subroutine MUMPS_TYPEINSSARBR(inode)
      implicit none
      integer, intent(in)::inode
      integer in
      cv_nodetype(inode)=-1
      in=cv_fils(inode)
      do while (in>0)
         in=cv_fils(in)
      end do
      in=-in
      do while(in.gt.0)
         call MUMPS_TYPEINSSARBR(in)
         in=cv_frere(in)
      enddo
      end subroutine MUMPS_TYPEINSSARBR
      subroutine MUMPS_WORKMEM_IMBALANCE(workload,memused,
     &                                    maxwork,minwork,maxmem,minmem)
      implicit none
      DOUBLE PRECISION,dimension(:),intent(in)::workload,
     &                                                  memused
      DOUBLE PRECISION,intent(out)::maxwork,minwork,maxmem,minmem
      intrinsic maxval,minval
      maxwork=maxval(workload)
      minwork=minval(workload, mask= workload > dble(0))
      maxmem=maxval(memused)
      minmem=minval(memused, mask= memused > dble(0))
      end subroutine MUMPS_WORKMEM_IMBALANCE
      subroutine MUMPS_FIX_ACCEPTED_MASTER(layernumber,nodenumber)
      implicit none
      integer layernumber,nodenumber
      integer i
      integer inode
      integer current_max,current_proc
      current_max = 0
      score = 0
      allowed_nodes = .FALSE.
      inode=cv_layer_p2node(layernumber)%t2_nodenumbers(nodenumber)
      do i=1,cv_layer_p2node(layernumber)%t2_cand(nodenumber,
     &     cv_slavef+1)
         current_proc=cv_layer_p2node(layernumber)%t2_cand(nodenumber,i)
         if ( current_proc .ge. 0) then
            score(mem_distribmpi(current_proc)) =
     &           score(mem_distribmpi(current_proc)) + 1
         endif
      enddo
      current_proc = cv_procnode(inode) - 1
      score(mem_distribmpi(current_proc)) =
     &     score(mem_distribmpi(current_proc)) + 1
      do i=0,nb_arch_nodes - 1
         if ( score(i) .gt. current_max ) then
            current_max = score(i)
            allowed_nodes = .FALSE.
            allowed_nodes(i) = .TRUE.
         else
            if(score(i) .eq. current_max) then
               allowed_nodes(i) = .TRUE.
            endif
         endif
      enddo
      return
      end subroutine MUMPS_FIX_ACCEPTED_MASTER
      end subroutine MUMPS_DISTRIBUTE
      subroutine MUMPS_RETURN_CANDIDATES(par2_nodes,cand,
     &                istat)
      integer, intent(out) :: par2_nodes(cv_nb_niv2), istat
      integer, intent(out) :: cand(:,:)
      character (len=48):: subname
      integer iloop
      istat=-1
      subname='MUMPS_RETURN_CANDIDATES'
      par2_nodes=cv_par2_nodes
      do iloop=1, cv_slavef+1
        cand(iloop,:)=cv_cand(:,iloop)
      enddo
      deallocate(cv_par2_nodes,cv_cand,STAT=istat)
      if(istat.ne.0) then
         if(cv_lp.gt.0)
     &   write(cv_lp,*)'Memory deallocation error in ',subname
         istat = cv_error_memdeloc
         return
      end if
      istat = 0
      return
      end subroutine MUMPS_RETURN_CANDIDATES
      subroutine MUMPS_INIT_ARCH_PARAMETERS(
     &     total_comm,working_comm,keep69,par,
     &     nbslaves,mem_distrib,informerr)
      implicit none
      include 'mpif.h'
      integer nbslaves
      integer, dimension(0:) :: mem_distrib
      integer total_comm,working_comm,keep69,par
      integer, dimension(:) ::informerr 
      integer myrank
      integer host,i,ierr
      integer,dimension(:),allocatable :: buffer_memdistrib
      ierr = 0
      myrank = -1
      host = -1
      ke69 = keep69
      cv_slavef = nbslaves
      if (ke69 .eq. 1) then
         return
      endif
      if ( allocated(mem_distribtmp) ) deallocate(mem_distribtmp )
      allocate( mem_distribtmp( 0:cv_slavef-1 ),
     &          buffer_memdistrib( 0:cv_slavef-1 ), stat=ierr )
      if ( ierr .gt. 0 ) then
         if(cv_mp.gt.0) write(cv_mp,*) 'pb allocation mem_dist'
         informerr(1) = -13
         informerr(2) = cv_slavef
         return
      end if
      mem_distribtmp = -1
      call MPI_COMM_RANK( total_comm, host, ierr )
      if ((par .eq. 1) .or. (host .ne. 0)) then
         call MPI_COMM_RANK( working_comm, myrank, ierr )
         call MUMPS_COMPUTE_DISTRIB(ierr,myrank,
     &        working_comm,mem_distrib)
         if ( ierr .ne. 0 ) then
            if(cv_mp.gt.0)
     &      write(cv_mp,*) 'pb in mumps_init_arch_parameters'
            informerr(1) = -13
            informerr(2) = cv_slavef
            return
         end if
         mem_distribtmp = mem_distrib
         call MUMPS_FIX_NODE_MASTER(ierr)
         if ( ierr .ne. 0 ) then
            if(cv_mp.gt.0) write(cv_mp,*)
     &'pb in mumps_init_arch_parameters'
            informerr(1) = -13
            informerr(2) = cv_slavef
            return
         endif
      endif
      if(ke69 .le. 0) then
         deallocate(mem_distribtmp)
         deallocate(buffer_memdistrib)
         return
      endif
      call MPI_ALLREDUCE(mem_distribtmp(0),buffer_memdistrib(0),
     &     cv_slavef,MPI_INTEGER,
     &     MPI_MAX,total_comm,ierr)
      mem_distribtmp = buffer_memdistrib
      deallocate (buffer_memdistrib)
      call MUMPS_COMPUTE_NB_ARCH_NODES()
        if((cv_slavef/nb_arch_nodes) .le. 4) then
          do i = 0, cv_slavef-1
            if ( mem_distrib(i) .NE. 1 ) then
              mem_distrib(i)=max(ke69/2,2)
            endif
          enddo
        endif
        if((nb_arch_nodes .eq. 1) .or.
     &     (nb_arch_nodes .eq. cv_slavef)) then
         ke69 = 1
         keep69 = 1
         deallocate(mem_distribtmp)
         return
      endif
      if (host .eq. 0) then
         if ( allocated(mem_distribmpi) ) deallocate(mem_distribmpi )
         allocate( mem_distribmpi( 0:cv_slavef-1 ), stat=ierr )
         if ( ierr .gt. 0 ) then
            if(cv_mp.gt.0) write(cv_mp,*) 'pb allocation mem_dist'
            informerr(1) = -13
            informerr(2) = cv_slavef
            return
         endif
         call MUMPS_ALLOC_ALLOW_MASTER(ierr)
         if(ierr .ne. 0 ) then
            return
         endif
         mem_distribmpi = mem_distribtmp
         call MUMPS_FIX_TABLE_OF_PROCESS(ierr)
         if ( ierr .ne. 0 ) then
            if(cv_mp.gt.0)
     &      write(cv_mp,*) 'pb in mumps_init_arch_parameters'
            informerr(1) = -13
            informerr(2) = cv_slavef
            return
         endif
      else
         deallocate(mem_distribtmp)
      endif
      return
      end subroutine MUMPS_INIT_ARCH_PARAMETERS
      subroutine MUMPS_COMPUTE_NB_ARCH_NODES()
      implicit none
      integer i
      nb_arch_nodes = 0
      do i=0,cv_slavef-1
         if(mem_distribtmp(i) .eq. i) then
            nb_arch_nodes = nb_arch_nodes + 1
         endif
      enddo
      return
      end subroutine MUMPS_COMPUTE_NB_ARCH_NODES
      subroutine MUMPS_FIX_TABLE_OF_PROCESS(ierr)
      implicit none
      external MUMPS_SORT_INT
      integer i,precnode,nodecount
      integer sizesmp
      integer ierr
      ierr = 0
      sizesmp = 0
      if ( allocated(table_of_process) )
     &     deallocate(table_of_process  )
      allocate( table_of_process(0:cv_slavef-1), stat=ierr )
      if ( ierr .gt. 0 ) then
         if(cv_mp.gt.0) write(cv_mp,*)
     &   'pb allocation in MUMPS_FIX_TABLE_OF_PROCESS'
         return
      end if
      do i=0,cv_slavef - 1
         table_of_process(i) = i
      enddo
      call MUMPS_SORT_INT(cv_slavef,mem_distribtmp(0),
     &     table_of_process(0))
      precnode = 0
      nodecount = 0
      do i=0,cv_slavef-1
         if(mem_distribtmp(i) .eq. precnode) then
            sizesmp = sizesmp + 1
            mem_distribtmp(i) = nodecount
            mem_distribmpi(table_of_process(i)) = nodecount
         else
            score(nodecount) = sizesmp
            sizesmp = 1
            nodecount = nodecount + 1
            precnode = mem_distribtmp(i)
            mem_distribtmp(i) = nodecount
            mem_distribmpi(table_of_process(i)) = nodecount
         endif
      enddo
      score(nodecount) = sizesmp
      do i=0,cv_slavef-1
         mem_distribtmp(i) = score(mem_distribtmp(i))
      enddo
      CALL MUMPS_SORT_INT_DEC(cv_slavef,mem_distribtmp(0),
     &     table_of_process(0))
      ierr = 0
      return
      end subroutine MUMPS_FIX_TABLE_OF_PROCESS
      subroutine MUMPS_FIX_NODE_MASTER(ierr)
      implicit none
      integer i,j,ierr
      integer idmaster
      idmaster = -1
      ierr = 0
      do i=0,cv_slavef-1
         if (mem_distribtmp(i) .eq. 1) then
            idmaster = i
            do j=i,cv_slavef-1
               if (mem_distribtmp(j) .eq. 1) then
                  mem_distribtmp(j) = idmaster
               else
                  mem_distribtmp(j) = 0
               endif
            enddo
            return
         else
            mem_distribtmp(i) = 0
         endif
      enddo
      if(cv_mp.gt.0) write(cv_mp,*)'problem in MUMPS_FIX_NODE_MASTER:
     &     cannot find a master'
      ierr = 1
      return
      end subroutine MUMPS_FIX_NODE_MASTER
      subroutine MUMPS_COMPUTE_DISTRIB(ierr,myrank,working_comm,
     &     mem_distrib)
      implicit none
      include 'mpif.h'
      integer ierr,resultlen,myrank,i,working_comm
      integer , dimension(0:) :: mem_distrib
      integer allocok
      character(len=MPI_MAX_PROCESSOR_NAME) name
      integer, dimension(:),allocatable :: namercv
      integer, dimension(:),allocatable :: myname
      integer lenrcv
      external MUMPS_COMPARE_TAB
      logical MUMPS_COMPARE_TAB
      ierr = 0
      call MPI_GET_PROCESSOR_NAME(name,resultlen,ierr)
      allocate(myname(resultlen),stat=allocok)
      if ( allocok .gt. 0 ) then
        if(cv_mp.gt.0) write(cv_mp,*)
     &  'pb allocation in compute_dist for myname'
        ierr = 1
        return
      end if
      do i=1, resultlen
         myname(i) = ichar(name(i:i))
      enddo
      do i=0, cv_slavef-1
         if(myrank .eq. i) then
            lenrcv = resultlen
         else
            lenrcv = 0
         endif
         call MPI_BCAST(lenrcv,1,MPI_INTEGER,i,
     &           working_comm,ierr)
         allocate(namercv(lenrcv),stat=allocok)
         if ( allocok .gt. 0 ) then
            if(cv_mp.gt.0) write(cv_mp,*)
     &      'pb allocation in compute_dist for namercv'
            ierr = 1
            return
         end if
         if(myrank .eq. i) then
            namercv = myname
         endif
         call MPI_BCAST(namercv,lenrcv,MPI_INTEGER,i,
     &        working_comm,ierr)
         if( MUMPS_COMPARE_TAB(myname,namercv,
     &        resultlen,lenrcv)) then
            mem_distrib(i)=1
         else
            mem_distrib(i)=ke69
         endif
         deallocate(namercv)
      enddo
      deallocate(myname)
      ierr = 0
      return
      end subroutine MUMPS_COMPUTE_DISTRIB
      subroutine MUMPS_GET_IDP1_PROC(current_proc,idarch,ierr)
      implicit none
      integer current_proc
      integer idarch,ierr
      ierr = 0
      if (current_proc .ge. cv_slavef) then
         ierr = -1
         return
      endif
      if (current_proc .lt. 0) then
         idarch = 1
         return
      else
         idarch = table_of_process(current_proc) + 1
      endif
      return
      end subroutine MUMPS_GET_IDP1_PROC
      subroutine MUMPS_END_ARCH_CV()
      if (allocated(table_of_process)) deallocate(table_of_process)
      if (allocated(allowed_nodes)) deallocate(allowed_nodes)
      if (allocated(score)) deallocate(score)
      if (allocated(mem_distribtmp)) deallocate(mem_distribtmp)
      if (allocated(mem_distribmpi)) deallocate(mem_distribmpi)
      return
      end subroutine MUMPS_END_ARCH_CV
      subroutine MUMPS_ALLOC_ALLOW_MASTER(ierr)
      integer ierr
      ierr = 0
      if (allocated(allowed_nodes)) deallocate(allowed_nodes)
      allocate( allowed_nodes(0:nb_arch_nodes-1),stat=ierr)
      if ( ierr .gt. 0 ) then
         if(cv_mp.gt.0) write(cv_mp,*)
     &   'pb allocation MUMPS_ALLOC_ALLOW_MASTER'
         ierr = -13
         return
      end if
      allowed_nodes = .FALSE.
       if (allocated(score)) deallocate(score)
      allocate( score(0:nb_arch_nodes-1),stat=ierr)
      if ( ierr .gt. 0 ) then
         if(cv_mp.gt.0) write(cv_mp,*)
     &   'pb allocation MUMPS_ALLOC_ALLOW_MASTER'
         ierr = -13
         return
      end if
      score = 0
      ierr = 0
      return
      end subroutine MUMPS_ALLOC_ALLOW_MASTER
      SUBROUTINE MUMPS_SORT_MMERGE(start1st,end1st,dim1,
     &                             start2nd,end2nd,dim2,
     &                             indx,
     &                             val, istat)
      implicit none
      integer, intent(in):: start1st,end1st,dim1,start2nd,end2nd,dim2
      integer, intent(inout):: indx(:) 
      DOUBLE PRECISION, intent(inout):: val(:) 
      INTEGER, intent(out) :: istat
      INTEGER, ALLOCATABLE, DIMENSION(:) :: index
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dummy1
      integer :: a,b,c
      integer :: allocok
      character (len=48):: subname
      subname = "MUMPS_SORT_MMERGE"
      istat=-1
      ALLOCATE(index(dim1+dim2),dummy1(dim1+dim2),stat=allocok)
      if ( allocok .gt. 0 ) then
         cv_info(1) = cv_error_memalloc
         cv_info(2) = dim1+dim2+dim1+dim2
         istat = cv_error_memalloc
         if(cv_lp.gt.0)
     &        write(cv_lp,*)
     &        'memory allocation error in ',subname
         return
      end if
      a=start1st
      b=start2nd
      c=1
      do while((a.LT.end1st+1).AND.(b.LT.end2nd+1))
         if(val(a).GT.val(b))then
            index(c)=indx(a)
            dummy1(c)=val(a)
            a=a+1
            c=c+1
         else
            index(c)=indx(b)
            dummy1(c)=val(b)
            b=b+1
            c=c+1
         endif
      end do
      if(a.LT.end1st+1) then
         do while(a.LT.end1st+1)
            index(c)=indx(a)
            dummy1(c)=val(a)
            a=a+1
            c=c+1
         enddo
      elseif(b.LT.end2nd+1) then
         do while(b.LT.end2nd+1)
            index(c)=indx(b)
            dummy1(c)=val(b)
            b=b+1
            c=c+1
         enddo
      endif
      indx(start1st:end1st)=index(1:dim1)
      val(start1st:end1st)=dummy1(1:dim1)
      indx(start2nd:end2nd)=index(dim1+1:dim1+dim2)
      val(start2nd:end2nd)=dummy1(dim1+1:dim1+dim2)
      DEALLOCATE(index,dummy1)
      istat=0
      return
      end SUBROUTINE MUMPS_SORT_MMERGE
      SUBROUTINE MUMPS_SORT_MSORT(istat,dim,indx,val1,val2)
      implicit none
      integer, intent(in):: dim
      integer, intent(inout):: indx(:)
      integer, intent(out)::istat
      DOUBLE PRECISION, intent(inout):: val1(:)
      DOUBLE PRECISION, intent(inout),optional:: val2(:)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: index, dummy1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dummy2
      integer, parameter :: ss = 35
      integer :: a,b,c,i,k,l,r,s,stackl(ss),stackr(ss)
      integer :: allocok
      character (len=48):: subname
      istat=-1
      subname = "MUMPS_SORT_MSORT"
      ALLOCATE(index(dim),dummy1(dim),dummy2(dim),stat=allocok)
      if (allocok.gt.0) then
        cv_info(1) = cv_error_memalloc
        cv_info(2) = 3*dim
        istat = cv_error_memalloc
        if(cv_lp.gt.0)
     &  write(cv_lp,*)'memory allocation error in ',subname
        return
      end if
      do i=1,dim
         index(i)=i
      enddo
      s = 1
      stackl(1) = 1
      stackr(1) = dim
 5511 CONTINUE
      l = stackl(s)
      r = stackr(s)
      k = (l+r) / 2
      if(l.LT.k) then
         if(s.GE.ss) stop 'maxsize of stack reached'
         s = s + 1
         stackl(s) = l
         stackr(s) = k
         goto 5511
      endif
 5512 CONTINUE
      l = stackl(s)
      r = stackr(s)
      k = (l+r) / 2
      if(k+1.LT.r) then
         if(s.GE.ss) stop 'maxsize of stack reached'
         s = s + 1
         stackl(s) = k+1
         stackr(s) = r
         goto 5511
      endif
 5513 CONTINUE
      l = stackl(s)
      r = stackr(s)
      k = (l+r) / 2
      a=l
      b=k+1
      c=1
      do while((a.LT.k+1).AND.(b.LT.r+1))
         if(val1(index(a)).GT.val1(index(b)))then
            dummy1(c)=index(a)
            a=a+1
            c=c+1
         else
            dummy1(c)=index(b)
            b=b+1
            c=c+1
         endif
      end do
      if(a.LT.k+1) then
         dummy1(c:r-l+1)=index(a:k)
      elseif(b.LT.r+1) then
         dummy1(c:r-l+1)=index(b:r)
      endif
      index(l:r)=dummy1(1:r-l+1)
      if(s.GT.1) then
         s = s - 1
         if(l.EQ.stackl(s)) goto 5512
         if(r.EQ.stackr(s)) goto 5513
      endif
      do i=1,dim
         dummy1(i)=indx(index(i))
      enddo
      indx=dummy1
      do i=1,dim
         dummy2(i)=val1(index(i))
      enddo
      val1=dummy2
      if(present(val2)) then
         do i=1,dim
            dummy2(i)=val2(index(i))
         enddo
         val2=dummy2
      endif
      istat=0
      DEALLOCATE(index,dummy1,dummy2)
      return
      end subroutine MUMPS_SORT_MSORT
      END MODULE MUMPS_STATIC_MAPPING
      SUBROUTINE MUMPS_SELECT_K38K20(N, SLAVEF, MP,
     &                   ICNTL13, KEEP, FRERE, ND, ISTAT)
      IMPLICIT NONE
      INTEGER, intent(in) :: N, SLAVEF, ICNTL13, MP
      INTEGER KEEP(500)
      INTEGER FRERE(N), ND(N)
      INTEGER, intent(out) :: ISTAT
      INTEGER IROOTTREE, SIZEROOT, NFRONT, I
      ISTAT = 0
      IF (KEEP(60).EQ.2 .or. KEEP(60).EQ.3 ) THEN
      ELSE
        IF((SLAVEF.EQ.1).OR.(ICNTL13.GT.0).OR.
     &     (KEEP(60).NE.0)) THEN
          KEEP(38) = 0
        ELSE
         IROOTTREE=-1
         SIZEROOT=-1
         DO I=1,N
            IF (FRERE(I).EQ.0) THEN
               NFRONT = ND(I)
               IF (NFRONT .GT.SIZEROOT) THEN
                  IROOTTREE = I
                  SIZEROOT  = NFRONT
               END IF
            END IF
         END DO
         IF ((IROOTTREE.EQ.-1).OR.(SIZEROOT.EQ.-1)) THEN
            ISTAT = -1
            RETURN
         ENDIF
         IF (SIZEROOT.LE.SLAVEF) THEN
            KEEP(38) = 0
         ELSE IF((SIZEROOT.GT.KEEP(37))
     &           .AND. (KEEP(53).EQ.0)
     &           ) THEN
            IF (MP.GT.0) WRITE(MP,*) 'A root of estimated size ',
     &           SIZEROOT,' has been selected for Scalapack.'
            KEEP(38) = IROOTTREE
         ELSE
            KEEP(38) = 0
             IF (MP.GT.0) WRITE(MP,'(A,I9,A)')
     &          ' WARNING: Largest root node of size ', SIZEROOT,
     &          ' not selected for parallel execution'
         END IF
         IF ((KEEP(38).EQ.0).AND.(KEEP(53).NE.0)) THEN
            KEEP(20) = IROOTTREE
         ELSE IF (KEEP(60).EQ.0) THEN
            KEEP(20) = 0
         ENDIF
       ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_SELECT_K38K20
      SUBROUTINE MUMPS_SPLITNODE_INTREE(inode,nfront,npiv,k,
     &     lnpivsplit, npivsplit, keep, n, fils, frere,
     &     nfsiz, ne, info5_nfrmax, k28_nsteps, nodetype,
     &     istat
     &                      , SIZEOFBLOCKS, LSIZEOFBLOCKS
     &                      , BLKON
     &      )
      implicit none
      integer, intent(in)::nfront,npiv
      integer, intent(in):: k
      integer, intent(in)::lnpivsplit
      integer, intent(in)::npivsplit(lnpivsplit)
      integer, intent(in):: inode
      integer, intent(out)::istat
      integer, intent(inout):: keep(500)
      integer, intent(inout):: k28_nsteps
      integer, intent(in)  :: info5_nfrmax
      integer, intent(in)  :: n
      integer, intent(inout)::frere(n), fils(n), nfsiz(n), ne(n)
      integer, intent(inout):: nodetype(n)
      integer, intent(in) :: LSIZEOFBLOCKS
      integer, intent(in) :: SIZEOFBLOCKS(LSIZEOFBLOCKS)
      logical,intent(in)  :: BLKON
      integer i,lev,in,in_son,in_father,in_grandpa,npiv_father,
     &     npiv_son,nfrontk,npivk,d1,f1,e1,dk,fk,next_father
      integer::ison,ifather
      character (len=48):: subname
      integer, parameter:: tsplit_beg=4 
      integer, parameter:: tsplit_mid=5 
      integer, parameter:: tsplit_last=6 
      istat=-1
      subname='SPLITNODE_INTREE'
      ison=-1
      ifather=-1
      nfrontk = nfront
      npivk = npiv
      npiv_son = npivsplit(1)
      keep(2)=max(keep(2),nfront-npiv_son)
      d1 = inode
      f1 = d1
      e1 = frere(d1)
      if (BLKON) then
         i= SIZEOFBLOCKS(f1)
         do while (i.lt.npiv_son)
            f1 = fils(f1)
             i = i + SIZEOFBLOCKS(f1)
         enddo
      else
         do i=1,npiv_son-1
           f1 = fils(f1)
         enddo
       endif
      ison = d1
      in_son = f1
      next_father = fils(in_son)
      do lev = 1, k-1
         ifather = next_father
         in_father = ifather
         npiv_son=   abs(npivsplit(lev))
         npiv_father=abs(npivsplit(lev+1))
         if (BLKON) then
            i= SIZEOFBLOCKS(in_father)
            do while (i.lt.npiv_father)
             in_father=fils(in_father)
             i = i + SIZEOFBLOCKS(in_father)
            enddo
         else
            do i=1,npiv_father-1
               in_father=fils(in_father)
            enddo
         endif
         frere(ison)=-ifather
         next_father = fils(in_father)
         fils(in_father)=-ison
         nfsiz(ison)=nfrontk  
         nfsiz(ifather)=nfrontk-npiv_son
         ne(ifather)=1
         keep(61)=keep(61)+1 
         IF (keep(79).EQ.0) THEN
           if( nfront-npiv_son > keep(9)) then
              nodetype(ifather) = 2
           else
              nodetype(ifather) = 1
           endif
         ELSE
          if (lev.EQ.1) then
             nodetype(ison) = tsplit_beg
          endif
          if (lev.eq.k-1) then
              nodetype(ifather) = tsplit_last  
          else
              nodetype(ifather) = tsplit_mid
          endif
          if  (npivsplit(lev+1) < 0) then
            if (lev.eq.k-1) then
               nodetype(ifather)=-tsplit_last
            else
              nodetype(ifather)=-tsplit_mid
            endif
          endif
         ENDIF
         nfrontk = nfrontk-npiv_son
         npivk = npivk - npiv_son
         ison = ifather
         in_son = in_father
      enddo
      dk = ifather
      fk = in_father
# if (check_mumps_static_mapping >= 3)
         write(6,*) ' Last (close to root) node in chain :', ifather
#endif 
      fils(f1) = next_father 
      frere(dk) = e1
      in = e1
      do while (in.gt.0)
         in=frere(in)
      end do
      in = -in
      do while(fils(in).gt.0)
         in=fils(in)
      end do
      in_grandpa = in
      if(fils(in_grandpa).eq.-d1) then
         fils(in_grandpa)=-dk
      else
         in=-fils(in_grandpa)
         do while(frere(in) .ne. d1)
            in=frere(in)
         end do
         frere(in) = dk
      end if
      k28_nsteps = k28_nsteps + k-1
      istat = 0
      return
      END SUBROUTINE MUMPS_SPLITNODE_INTREE
      subroutine MUMPS_SETUP_CAND_CHAIN(n, nb_niv2,
     &        frere, nodetype, par2_nodes,
     &        procnode, cand, inode_chain, slavef, dummy, nbcand, istat)
      implicit none
      integer, intent(in) :: n, nb_niv2, slavef
      integer,intent(in)::frere(n)
      integer, intent(inout) :: par2_nodes(nb_niv2), procnode(n)
      integer,intent(inout)::nodetype(n)
      integer,intent(inout)::cand(nb_niv2, slavef+1)
      integer,intent(in)::inode_chain
      integer,intent(inout)::dummy, nbcand
      integer,intent(out):: istat
      integer, parameter:: tsplit_beg=4 
      integer, parameter:: tsplit_mid=5 
      integer, parameter:: tsplit_last=6 
      integer, parameter:: invalid=-9999
      integer            :: inode, ifather, k
      logical            :: last_iteration_reached
      istat = -1
      inode = inode_chain
      k  = 1  
      do
         if (.not. (frere(inode) .lt. 0) ) then
          write(*,*) " Internal error 0 in SETUP_CAND",
     &    frere(inode), inode
          CALL MUMPS_ABORT()
        endif
        ifather = -frere(inode)
        last_iteration_reached = (abs(nodetype(ifather)).eq.tsplit_last)
        par2_nodes(dummy+1) = ifather
        procnode(ifather) = cand(dummy,1) + 1
        if  ( (nodetype(ifather).eq.tsplit_mid) .or.
     &        (nodetype(ifather).eq.tsplit_last) ) then
           if (nbcand.lt.2) then
            par2_nodes(dummy+1) = ifather
            procnode(ifather) = procnode(inode) 
            cand(dummy+1,:) = cand(dummy,:)
            dummy = dummy + 1
            write(6,*) ' Mapping property', 
     &                 ' of procs in chain lost '
            CALL MUMPS_ABORT()
          endif
          cand(dummy+1,1:nbcand-1+k-1) = cand(dummy,2:nbcand+k-1)
          cand(dummy+1,nbcand-1+k) = procnode(inode)-1 
          cand(dummy+1,nbcand-1+k+1:slavef) = invalid
          nbcand = nbcand -1
          k = k + 1
        else if ( (nodetype(ifather).eq.-tsplit_mid) .or.
     &            (nodetype(ifather).eq.-tsplit_last) ) then
          if (nodetype(inode).eq.tsplit_beg) then
            nodetype(inode)=2
          else
            nodetype(inode)=tsplit_last
          endif
          if (nodetype(ifather) .eq. -tsplit_last) then
            nodetype(ifather) = 2 
          else
            nodetype(ifather) = tsplit_beg
          endif
          cand(dummy+1,1:nbcand-1+k-1) = cand(dummy,2:nbcand+k-1)
          cand(dummy+1,nbcand-1+k) = procnode(inode)-1 
          nbcand = nbcand+k-1
          k = 1
       else
          write(6,*) ' Internal error 2 in SETUP_CAND', 
     &    ' in, ifather =', inode, ifather,
     &    ' nodetype(ifather) ', nodetype(ifather)
          CALL MUMPS_ABORT()
        endif
        cand(dummy+1,slavef+1)= nbcand
        dummy = dummy+1
        if  (last_iteration_reached) exit
        inode = ifather 
      end do
      istat = 0
      end subroutine MUMPS_SETUP_CAND_CHAIN
      subroutine MUMPS_GET_SPLIT_4_PERF(inode, nfront, npiv, nproc,
     &                                  k, lnpivsplit, npivsplit,
     &                                  n, frere, keep,
     &                                  fils, BLKON, sizeofblocks,
     &                                  istat)
      implicit none
      integer,intent(in)::inode, nfront, npiv, lnpivsplit, n
      integer,intent(in)::frere(n)
      integer,intent(in) :: fils(n)
      logical, intent(in) :: BLKON
      integer, intent(in) :: sizeofblocks(*)
      integer,intent(in)::keep(500)
      double precision, intent(in):: nproc
      integer,intent(out)::k, npivsplit(lnpivsplit), istat
      logical :: nosplit
      integer :: inode_tmp
      integer :: kk, optimization_strategy, nass, npiv2
      double precision :: nproc2
      integer :: npivOld, npivNew
      double precision :: timeFacOld, timeFacNew, timeAss
      double precision ,parameter :: alpha=8.0D9
      double precision ,parameter :: gamma=1.2D9
      nosplit = npiv .le. npiv4equilibreRows(nfront, nproc)
      optimization_strategy = 0 
      nosplit = nosplit .or. (frere(inode) .eq. 0)
      if ( nosplit ) then
         k = 1
         npivsplit(1) = npiv
         istat = 0
         return
      endif
      if (nproc .le. 1.0d0) then 
         k = 1
         npivsplit(1) = npiv
         istat = -1
         return
      endif
      nproc2 = nproc
      nass = 0
      kk = 0
      inode_tmp = inode
      do while (nass .lt. npiv)
        if ((nproc2 .eq. 2.0d0) .or.
     &           (nfront - nass .le. 6*keep(9))) then
          npiv2 = npiv - nass
        else if (nproc2 .gt. 2) then
          if (optimization_strategy .eq. 0) then 
            npiv2 = min(npiv - nass,
     &                  npiv4equilibreRows(nfront - nass, nproc2 ))
          else if (optimization_strategy .eq. 1) then 
            if (nproc2 .eq. nproc) then
              npiv2 = min(npiv - nass,
     &                    npiv4equilibreFlops(nfront - nass, nproc2 ))
            else 
              npiv2 = min(npiv - nass,
     &                    npiv4equilibreRows(nfront - nass, nproc2 ))
            endif
          else
            write(*,*) "Internal error in MUMPS_GET_SPLIT_4_PERF,"
            write(*,*) "optimization_strategy not implemented"
            call MUMPS_ABORT()
          endif
        endif
        kk = kk + 1
        IF (BLKON) THEN
          npivsplit(kk) = 0
          DO WHILE (npivsplit(kk) .LT. npiv2 .and. inode_tmp .gt. 0)
            npivsplit(kk) = npivsplit(kk) + sizeofblocks(inode_tmp)
            inode_tmp= fils(inode_tmp)
          ENDDO
          npiv2 = npivsplit(kk)
        ELSE
          npivsplit(kk) = npiv2
        ENDIF
        if (keep(79) .ge. 1       
     &    .and. kk .ne. 1) then   
          if (optimization_strategy .eq. 0) then 
            npivOld = min(npiv - nass,
     &              npiv4equilibreRows(nfront - nass, nproc ))
            npivNew = min(npiv - nass,
     &              npiv4equilibreRows(nfront - nass, nproc2 - 1.0d0))
          else if (optimization_strategy .eq. 1) then 
            npivOld = min(npiv - nass,
     &              npiv4equilibreFlops(nfront - nass, nproc ))
            npivNew = min(npiv - nass,
     &              npiv4equilibreRows(nfront - nass, nproc2 - 1.0d0))
          else
            write(*,*) "Internal error in MUMPS_GET_SPLIT_4_PERF,"
            write(*,*) "optimization_strategy not implemented"
            call MUMPS_ABORT()
          endif
          timeAss = timeAssembly(int(nfront-nass,8), nproc2)
          timeFacOld = timeFacto(int(nfront-nass,8), int(npivOld,8),
     &                           nproc)
          timeFacNew = timeFacto(int(nfront-nass,8),int(npivNew,8),
     &                           nproc2-1)
          if ( (flopsFactoPanel(int(npivOld,8),int(nfront-nass,8))+
     &          flopsUpdate(int(nfront-nass-npivOld,8),
     &          int(nfront-nass-npivOld,8), int(npivOld,8)))/
     &          (timeFacOld+timeAss)
     &       .gt. (flopsFactoPanel(int(npivNew,8),int(nfront-nass,8))+
     &          flopsUpdate(int(nfront-nass-npivNew,8),
     &          int(nfront-nass-npivNew,8), int(npivNew,8)))/
     &          timeFacNew ) then 
            npivsplit(kk) = -npiv2
            nproc2 = nproc
          else
            nproc2 = nproc2 - 1.0d0
            npiv2 = npivNew
            npivsplit(kk)=npivNew
          endif
        endif
        nass = nass + npiv2
      enddo
      k = kk
      istat=0
      return
      CONTAINS   
      function npiv4equilibreRows(nfront, nproc)
      implicit none
      integer npiv4equilibreRows
      integer, intent(in) :: nfront
      double precision, intent(in) :: nproc
        npiv4equilibreRows = max(1, int(dble(nfront)/nproc))
        return 
      end function npiv4equilibreRows
      function  npiv4equilibreFlops(nfront, nproc)
      implicit none
      integer npiv4equilibreFlops
      integer, intent(in) :: nfront
      double precision, intent(in) :: nproc
        double precision::n,s,a,b,c,sdelta,npiv
        n = dble(nfront)
        s = nproc - 1.0d0
        a = s/3.+1.
        b = -3.*n - s*n - s/2.
        c = 2.*n**2 + s*n + s/6.
        sdelta = (b*b) - 4*a*c
        if (sdelta < 0.0E0) then
          WRITE(*,*) "Delta < 0 in npiv4equilibreFlops"
          call MUMPS_ABORT()
        endif
        sdelta = sqrt(sdelta)
        npiv = (-b - sdelta)/(2*a)
        npiv4equilibreFlops = max(1, int(npiv))
        return
      end function npiv4equilibreFlops
      function flopsFactoPanel(nbrows, nbcols)
        integer(8) :: nbrows, nbcols
        double precision :: flopsFactoPanel
        flopsFactoPanel = (nbrows*((-1.d0/3.d0)*nbrows**2 +
     &          (nbcols + 1.d0/2.d0)*nbrows +
     &          (nbcols + 1.d0/6.d0)))
      end function flopsFactoPanel
      function flopsUpdate(m, n, k)
        integer(8) :: m, n, k
        double precision :: flopsUpdate
        flopsUpdate = dble(2*m*n*k + m*k**2)
      end function flopsUpdate
      function timeFacto(nfront, npiv, nproc)
        integer(8), intent(in) :: nfront, npiv
        double precision, intent(in) :: nproc
        double precision :: timeFacto
        timeFacto = (max(flopsFactoPanel(npiv,nfront),
     &              flopsUpdate(nfront-npiv, nfront-npiv, npiv)/
     &              (nproc-1))/alpha)
      end function timeFacto
      function timeNIV1(nfront, npiv)
        integer(8) :: nfront, npiv
        double precision :: timeNIV1
        timeNIV1 = ((flopsFactoPanel(npiv, nfront) +
     &           flopsUpdate(nfront - npiv, nfront - npiv, npiv))/alpha)
      end function timeNIV1
      function timeAssembly(n, p)
        integer(8) :: n
        double precision, intent(in) :: p
        double precision :: timeAssembly
        timeAssembly = ((n*n/p)/(gamma/(log(p)/log(2.0d0))))
      end function timeAssembly
      end subroutine MUMPS_GET_SPLIT_4_PERF
