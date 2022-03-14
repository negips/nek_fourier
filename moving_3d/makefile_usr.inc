
SRC_INC     := ./src_inc
SRC_USR     := ./src_usr
FRAME_SRC   := ./frame
SRC_ARN     := ./src_arn

$(OBJDIR)/frame.o               :$(FRAME_SRC)/frame.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrlog_block.o       :$(FRAME_SRC)/mntrlog_block.f;          $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrlog.o             :$(FRAME_SRC)/mntrlog.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrtmr_block.o       :$(FRAME_SRC)/mntrtmr_block.f;          $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrtmr.o             :$(FRAME_SRC)/mntrtmr.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/rprm_block.o          :$(FRAME_SRC)/rprm_block.f;             $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/rprm.o                :$(FRAME_SRC)/rprm.f;                   $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/io_tools.o            :$(FRAME_SRC)/io_tools.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/io_tools_block.o      :$(FRAME_SRC)/io_tools_block.f;         $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/math_tools.o          :$(FRAME_SRC)/math_tools.f;             $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/solve_f3d.o           :$(SRC_USR)/solve_f3d.f ;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/solvep_f3d.o          :$(SRC_USR)/solvep_f3d.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/weak_lapl_cyl.o       :$(SRC_USR)/weak_lapl_cyl.f;            $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/solver_f3d.o          :$(SRC_USR)/solver_f3d.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/frame_f3d.o           :$(SRC_USR)/frame_f3d.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/tstepper.o            :$(SRC_ARN)/tstepper.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/conht_tools.o         :$(SRC_ARN)/conht_tools.f;              $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_arp.o             :$(SRC_ARN)/arn_arp.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_arp_io.o          :$(SRC_ARN)/arn_arp_io.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/chkpoint.o            :$(SRC_ARN)/chkpoint.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/chkptd.o              :$(SRC_ARN)/chkptd.f;                   $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_cmath.o           :$(SRC_ARN)/arn_cmath.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fs_ale.o              :$(SRC_USR)/fs_ale.f;                   $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fs_global_surf.o      :$(SRC_USR)/fs_global_surf.f;           $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/frame_fs.o            :$(SRC_USR)/frame_fs.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/localfilter.o         :$(SRC_USR)/localfilter.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@



