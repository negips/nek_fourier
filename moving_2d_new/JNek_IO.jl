#     Port for reader_par.f
#     Author:     Prabal Negi
#

      module JNek_IO

      using MPI

      export read_re2_hdr,
             read_re2, 
             read_ma2

      function __init__()

        if MPI.Initialized() == false

          MPI.Init()

        end    
          
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)

        if rank == 0
          println("Initialied MPI in Module JNek_IO")
        end  

        return nothing
      end 

#----------------------------------------------------------------------

      function init()

#       Can't reinitialize MPI        

        return nothing
      end 

#----------------------------------------------------------------------


      function read_re2_hdr(fid::IOStream, rank)

        println("Re2: Reading Header on rank $(rank)")
        
        nbytes  = 20*4
        hdrutf  = read(fid,nbytes)
        hdrC    = Char.(hdrutf)
        version = String(hdrC[1:5])
        nelgtS  = String(hdrC[6:14])
        ldimrS  = String(hdrC[15:17])
        nelgvS  = String(hdrC[18:26])

        nelgt   = parse(Int, nelgtS)
        ldimr   = parse(Int, ldimrS)
        nelgv   = parse(Int, nelgvS)

        hdr     = String(hdrC)

        test    = read(fid,Float32)

        if_byte_swap = byte_swap_test(test)

        return hdr,version,nelgt,ldimr,nelgv,if_byte_swap
      end     # read_re2_hdr

#---------------------------------------------------------------------- 

      function read_re2(f::String, nid0::Int64)

        hdr = repeat(" ", 26)
        version = repeat(" ", 5)
        nelgt = 0
        nelgv = 0
        ldimr = 0

        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)

        if rank == nid0
          println("Reading $(f) on rank $(MPI.Comm_rank(comm))\n")

          fid = open(f, "r")
          hdr,version,nelgt,ldimr,nelgv,if_byte_swap = read_re2_hdr(fid,rank)
        end  

        MPI.Barrier(comm)

        buf = MPI.Buffer(hdr,26,MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        buf = MPI.Buffer(version,5,MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        nelgt = MPI.bcast(nelgt,     nid0,comm)
        ldimr = MPI.bcast(ldimr,     nid0,comm)
        nelgv = MPI.bcast(nelgv,     nid0,comm)

        wdsizi = 4
        if cmp(version,"#v002") == 0 || cmp(version, "#v003")
          wdsizi = 8
        end   

#       Read the Mesh data here
        xc,yc = read_re2_mesh(fid,nid0,ldimr,nelgt,wdsizi)

        ncurve,curveieg,curveiside,curveparam,curvetype = read_re2_curve(fid,nid0,ldimr,nelgt,wdsizi)

        cbl,bl = read_re2_bc(fid,nid0,ldimr,nelgt,wdsizi)

        if rank == nid0
          close(fid)
        end 

        return hdr,version,nelgt,ldimr,nelgv,xc,yc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl
      end     # read_re2

#---------------------------------------------------------------------- 

      function read_re2_mesh(fid::IOStream, nid0::Int64,ldim::Int64,nelgt::Int64,wdsizi::Int64)

#       Pointer to re2 data in file.
#       Header + test pattern byte length
        recpos  = 84
        seek(fid,recpos)

        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)

        nc   = 2^ldim
        xc   = Array{Float64,2}(undef,nc,nelgt)
        yc   = Array{Float64,2}(undef,nc,nelgt)

        len     = 1 + ldim*(2^ldim)             # group + 2x4 for 2d, 3x8 for 3d
        nbytes  = len*wdsizi                    
       
        if (wdsizi == 4)
          tmp  = Vector{Float32}(undef,len)
        else
          tmp  = Vector{Float64}(undef,len)
        end
        tmp64  = Vector{Float64}(undef,len)

        for e = 1:nelgt
          if rank == nid0
            read!(fid,tmp)
            tmp64 = tmp
          end        
          if ldim == 2
            group   = tmp64[1]
#           In principle I should broadcast the data to different processors 
            xc[1:4,e] = tmp64[2:5]
            yc[1:4,e] = tmp64[6:9]
          end
        end  

        return xc,yc
      end     # read_re2_mesh

#---------------------------------------------------------------------- 

      function read_re2_curve(fid::IOStream, nid0::Int64,ldim::Int64,nelt::Int64,wdsizi::Int64)

        curveieg   = Vector{Int64}(undef,1)
        curveiside = Vector{Int64}(undef,1)
        curveparam = Matrix{Float64}(undef,5,1)
        curvetype  = Vector{Char}(undef,1)
        ncurve::Int64 = 0

#       Pointer to re2 data in file.
#       Header + test pattern
#        recpos  = 84
#        if wdsizi == 4
#          meshlength = (1 + ldim*(2^ldim))*4
#        else
#          meshlength = (1 + ldim*(2^ldim))*8
#        end
#        recpos = recpos + meshlength*nelt 
#        seek(fid,recpos)

        if wdsizi == 4
          nc   = read(fid,Float32)
        else
          nc   = read(fid,Float64)
        end  


        ncurve = nc
        if ncurve > 0
          curveieg   = Vector{Int64}(undef,ncurve)
          curveiside = Vector{Int64}(undef,ncurve)
          curveparam = Matrix{Float64}(undef,5,ncurve)
          curvetype  = Vector{Char}(undef,ncurve)
        end 

#        len     = 2 + 1 + 5             # ieg iside curve(5) ccurve
        if (wdsizi == 4)
          tmpi  = Vector{Int32}(undef,2)
          tmpr  = Vector{Float32}(undef,5)
          tmpc  = Vector{Char}(undef,1)
        else
          tmpi  = Vector{Int64}(undef,2)
          tmpr  = Vector{Float64}(undef,5)
          tmpc  = Vector{Char}(undef,2)
        end
       
        for i in 1:nc
          
          read!(fid,tmpi)     # Read ieg, iside
          curveieg[i]   = tmpi[1]
          curveiside[i] = tmpi[2]
          read!(fid,tmpr)     # Read Curve params 
          curveparam[:,i] = tmpr            
          read!(fid,tmpc)     # Read Curve type (ccurve)
          curvetype[i]    = tmpc 
        end  

        return ncurve,curveieg,curveiside,curveparam,curvetype 
      end     # read_re2_curve

#---------------------------------------------------------------------- 
      function read_re2_bc(fid::IOStream, nid0::Int64,ldim::Int64,nelt::Int64,wdsizi::Int64)

        nbc::Int64 = 0

        nfaces = 2*ldim
        cbl = Array{String}(undef,nfaces,nelt)
        bl  = zeros(Float64,5,nfaces,nelt)

#       Initialize Array 
        for i in 1:nelt, j in 1:nfaces
          cbl[j,i] = "E  "
        end

        nbc   = read(fid,Float64)

        if (wdsizi == 4)
          tmpi  = Vector{Int32}(undef,2)
          tmpr  = Vector{Float32}(undef,5)
          tmpc  = Vector{Char}(undef,1)
        else
          tmpi  = Vector{Float64}(undef,2)
          tmpr  = Vector{Float64}(undef,5)
          tmpc  = Vector{Char}(undef,2)
        end

#       len     = 2 + 1 + 5        # eg iside bl(5) cbl
        e::Int64  = 0 
        f::Int64  = 0
      
        for i in 1:nbc
          read!(fid,tmpi)     # Read eg, iside
          e = tmpi[1]
          f = tmpi[2]
          read!(fid,tmpr)     # Read bl params 
          bl[:,f,e] = tmpr            
          read!(fid,tmpc)     # Read bl type ("E  ", etc)
          cbl[f,e]  = String(tmpc)
        end  

#       place holder
        return cbl,bl
      end     # read_re2_bc

#---------------------------------------------------------------------- 

      function read_ma2(f::String, nid0::Int64)

        comm = MPI.COMM_WORLD

        if MPI.Comm_rank(comm) == nid0
          println("Reading $(f) on rank $(MPI.Comm_rank(comm))")
        end
        
        MPI.Barrier(comm)

        rank = MPI.Comm_rank(comm)

        if rank == nid0
          fid = open(f, "r")
        end

        close(fid)

        return nothing
      end     # read_ma2

#----------------------------------------------------------------------

      function byte_swap_test(test::Float32)

        pattern::Float32  = 6.54321
        eps::Float32      = 0.00020 
        if_byte_swap      = false
         
        etest = abs(test - pattern)
        if (etest>eps) 
          if_byte_swap    = true
        end

        return if_byte_swap

      end  

#----------------------------------------------------------------------

      function read_fld(f::String, MPI::Module, nid0::Int64)


        comm = MPI.COMM_WORLD

        if MPI.Comm_rank(comm) == nid0
          println("Reading $(f) on rank $(MPI.Comm_rank(comm))\n")
        end
        
        MPI.Barrier(comm)

        fid = open(f, "r")

        rank = MPI.Comm_rank(comm)

        hdr,version,wdsizi,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,if_byte_swap = read_fld_std_hdr(fid,rank,nid0)

        buf = MPI.Buffer(hdr,length(hdr),MPI.CHAR)
        MPI.Bcast!(buf,       nid0,comm)

        buf = MPI.Buffer(version,length(version),MPI.CHAR)
        MPI.Bcast!(buf,       nid0,comm)

        wdsizi          = MPI.bcast(wdsizi,           nid0,comm)
        nx              = MPI.bcast(nx,               nid0,comm)
        ny              = MPI.bcast(ny,               nid0,comm)
        nz              = MPI.bcast(nz,               nid0,comm)
        nel             = MPI.bcast(nel,              nid0,comm)
        nelgt           = MPI.bcast(nelgt,            nid0,comm)
        time            = MPI.bcast(time,             nid0,comm)
        istep           = MPI.bcast(istep,            nid0,comm)
        fid0            = MPI.bcast(fid0,             nid0,comm)
        nfileo          = MPI.bcast(nfileo,           nid0,comm)

        buf             = MPI.Buffer(rdcode,length(rdcode),MPI.CHAR)
        MPI.Bcast!(buf,       nid0,comm)
        p0th            = MPI.bcast(p0th,             nid0,comm)
        ifprmesh        = MPI.bcast(ifprmesh,         nid0,comm)
        if_byte_swap    = MPI.bcast(if_byte_swap,     nid0,comm)

#       Read the data here        
        glnum,x,y,z,u,v,w,p,t = read_fld_data(fid, nid0,nx,ny,nz,nelgt,rdcode,wdsizi)

        close(fid)

        return hdr,version,wdsizi,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x,y,z,u,v,w,p,t
      end     # read_fld

#---------------------------------------------------------------------- 

      function read_fld_std_hdr(fid::IOStream, rank, nid0)

#        comm = MPI.COMM_WORLD
        
        hdr           = repeat(" ", 132)
        version       = repeat(" ", 5)
        wdsize        = 0 
        nx            = 0
        ny            = 0
        nz            = 0
        nel           = 0
        nelgt         = 0
        time          = 0.0
        istep         = 0
        fid0          = 0
        nfileo        = 0
        rdcode        = repeat(" ",10)
        p0th          = 0.0 
        ifprmesh      = false

        if rank == nid0
          println("Fld: Reading Header on rank $(rank)")
        
          nbytes        = 132
          hdrutf        = read(fid,nbytes)
          hdrC          = Char.(hdrutf)
          
          st            = 2         # step
          i             = 1
          j             = 4
          version       = String(hdrC[i:j])
          i             = j+st # 7 
          j             = i+0
          wdsizS        = String(hdrC[i:j])
          i             = j+st # 
          j             = i+1
          nxS           = String(hdrC[i:j])
          i             = j+st # 
          j             = i+1
          nyS           = String(hdrC[i:j])
          i             = j+st # 
          j             = i+1
          nzS           = String(hdrC[i:j])
          i             = j+st # 
          j             = i+9
          nelS          = String(hdrC[i:j])
          i             = j+st # 
          j             = i+9
          nelgS         = String(hdrC[i:j])
          i             = j+st # 
          j             = i+19
          timeS         = String(hdrC[i:j])
          i             = j+st # 
          j             = i+8
          istepS        = String(hdrC[i:j])
          i             = j+st # 
          j             = i+5
          fid0S         = String(hdrC[i:j])
          i             = j+st # 
          j             = i+5
          nfileoS       = String(hdrC[i:j])
          i             = j+st # 
          j             = i+9
          rdcodeS       = String(hdrC[i:j])
          i             = j+st # 
          j             = i+14
          p0thS         = String(hdrC[i:j])
          i             = j+st # 
          j             = i+0
          ifpr_meshS    = String(hdrC[i:j])
    
          wdsize        = parse(Int, wdsizS)
          nx            = parse(Int, nxS)
          ny            = parse(Int, nyS)
          nz            = parse(Int, nzS)
          nel           = parse(Int, nelS)
          nelgt         = parse(Int, nelgS)
          time          = parse(Float64, timeS)
          istep         = parse(Int, istepS)
          fid0          = parse(Int, fid0S)
          nfileo        = parse(Int, nfileoS)
          rdcode        = rdcodeS
          p0th          = parse(Float64, p0thS)
          if (ifpr_meshS == "T")
            ifprmesh = true
          else
            ifprmesh = false
          end  
          hdr           = String(hdrC) 
        end

        test    = read(fid,Float32)

        if_byte_swap = byte_swap_test(test)

        return hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,if_byte_swap
      end     # read_fld_std_hdr
#---------------------------------------------------------------------- 
      function read_fld_data(fid::IOStream, nid0::Int64,nx::Int64,ny::Int64,nz::Int64,nelgt::Int64,rdcode::String,wdsizi::Int64)

#       Pointer to re2 data in file.
#       Header + test pattern byte length
        recpos  = 132+4
        seek(fid,recpos)

        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)

        if (wdsizi==4)
          gnum     = Vector{Int32}(undef,nelgt)
        else
          gnum     = Vector{Int64}(undef,nelgt)
        end
        glnum      = Vector{Int64}(undef,nelgt)

        read!(fid,gnum)
        glnum = gnum

        ldim = 3
        if nz == 1
          ldim = 2
        end  

        nxyz    = nx*ny*nz
        len     = ldim*(nxyz)
        nbytes  = len*wdsizi                    
       
        if (wdsizi == 4)
          tmp   = Array{Float32}(undef,nx,ny,nz)
          tmpv  = Array{Float32}(undef,nx,ny,nz,ldim)
        else
          tmp   = Array{Float64}(undef,nx,ny,nz)
          tmpv  = Array{Float64}(undef,nx,ny,nz,ldim)
        end
        tmp64   = Array{Float64}(undef,nx,ny,nz)
        tmpv64  = Array{Float64}(undef,nx,ny,nz,ldim)

        ifxo      = false
        ifuo      = false
        ifpo      = false
        ifto      = false
        ifpso     = false
        nt        = 0
        nps       = 0

        i = 0
        for s in rdcode
          i = i+1
          if s=='X'
            ifxo = true
          end
          if s=='U'
            ifuo = true
          end
          if s=='P'
            ifpo = true
          end
          if s=='T'
            ifto = true
            nt   = 1
          end
          if s=='S'
            ifpso = true
            nps   = parse(Int64,rdcode[i+1:i+2])
          end
        end  

        nt = nt+nps
        x   = Array{Float64,4}(undef,nx,ny,nz,nelgt)
        y   = Array{Float64,4}(undef,nx,ny,nz,nelgt)
        z   = Array{Float64,4}(undef,nx,ny,nz,nelgt)

        u   = Array{Float64,4}(undef,nx,ny,nz,nelgt)
        v   = Array{Float64,4}(undef,nx,ny,nz,nelgt)
        w   = Array{Float64,4}(undef,nx,ny,nz,nelgt)

        p   = Array{Float64,4}(undef,nx,ny,nz,nelgt)

       
        t   = Array{Float64,5}(undef,nx,ny,nz,nelgt,nt)

        if (ifxo)
          for e = 1:nelgt
            if rank == nid0
              read!(fid,tmpv)
              tmpv64 = tmpv
            end
            gn = glnum[e]             
            if ldim == 2
#             In principle I should broadcast the data to different processors
              x[:,:,:,gn]  = tmpv64[:,:,:,1]
              y[:,:,:,gn]  = tmpv64[:,:,:,2]
            else
              x[:,:,:,gn]  = tmpv64[:,:,:,1]
              y[:,:,:,gn]  = tmpv64[:,:,:,2]
              z[:,:,:,gn]  = tmpv64[:,:,:,3]
            end
          end
        end  # ifxo  

        if (ifuo)
          for e = 1:nelgt
            if rank == nid0
              read!(fid,tmpv)
              tmpv64 = tmpv
            end        
            gn = glnum[e]             
            if ldim == 2
#             In principle I should broadcast the data to different processors
              u[:,:,:,gn]  = tmpv64[:,:,:,1]
              v[:,:,:,gn]  = tmpv64[:,:,:,2]
            else
              u[:,:,:,gn]  = tmpv64[:,:,:,1]
              v[:,:,:,gn]  = tmpv64[:,:,:,2]
              w[:,:,:,gn]  = tmpv64[:,:,:,3]
            end
          end
        end  # ifuo  

        if (ifpo)
          for e = 1:nelgt
            if rank == nid0
              read!(fid,tmp)
              tmp64 = tmp
            end        
#           In principle I should broadcast the data to different processors 
            gn = glnum[e]             
            p[:,:,:,gn]  = tmp64
          end
        end  # ifpo  

        if (ifto || ifpso)
          for i = 1:nt
            for e = 1:nelgt
              if rank == nid0
                read!(fid,tmp)
                tmp64 = tmp
              end        
#             In principle I should broadcast the data to different processors 
              gn = glnum[e]             
              t[:,:,:,e,i]  = tmp64
            end
          end 
        end  # ifto || ifpso 

       return glnum,x,y,z,u,v,w,p,t
      end     # read_fld_data

#---------------------------------------------------------------------- 

#---------------------------------------------------------------------- 
      end   # Module JNek_IO










