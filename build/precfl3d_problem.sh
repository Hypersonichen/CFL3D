#mpif90 -DP3D_SINGLE -DLINUX -DCGNS -DDBLE_PRECSN -I/home/ubuntu/cfl3d/cgns/include  -o precfl3d sizer.o get_cmd_args.o rpatch0.o outbuf.o getibk0.o global.o global2.o pointers.o compg2n.o umalloc.o cntsurf.o lead.o pre_patch.o pre_blockbc.o pre_embed.o pre_period.o usrint.o termn8.o getdhdr.o global0.o readkey.o parser.o ccomplex.o cgnstools.o setseg.o my_flush.o main.o  -L$(HOME)/cfl3d/cgns/lib -lcgns -L../../cfl/libs -lcommon
mpif90 -DP3D_SINGLE -DLINUX  -DDBLE_PRECSN   -o precfl3d sizer.o get_cmd_args.o rpatch0.o outbuf.o getibk0.o global.o global2.o pointers.o compg2n.o umalloc.o cntsurf.o lead.o pre_patch.o pre_blockbc.o pre_embed.o pre_period.o usrint.o termn8.o getdhdr.o global0.o readkey.o parser.o ccomplex.o cgnstools.o setseg.o my_flush.o main.o   -L../../cfl/libs -lcommon

