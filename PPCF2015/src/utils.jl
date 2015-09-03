using Grid

function resize_transfer_matrix{T<:FloatingPoint}(A::Array{T,2},x::Array{T,1},y::Array{T,1},nxx::Int,nyy::Int)

    nr, nc = size(A)

    nx = length(x)
    xrange = linrange(x[1],x[end],nx)
    newxrange = linrange(x[1],x[end],nxx)

    ny = length(y)
    yrange = linrange(y[1],y[end],ny)
    newyrange = linrange(y[1],y[end],nyy)

    Ap = zeros(nr,nxx*nyy)

    for i=1:nr
        d = reshape(A[i,:],ny,nx)
        di = Grid.CoordInterpGrid((yrange,xrange),d,Grid.BCnearest,Grid.InterpLinear);
        dp = zeros(nyy,nxx)

        for (xi, xp) in enumerate(newxrange)
            for (yi,yp) in enumerate(newyrange)
                dp[yi,xi] = di[yp,xp]
            end
        end

        Ap[i,:]=reshape(dp,1,int(nxx*nyy))
    end

    return Ap,[newxrange],[newyrange]
end

function
  resize_transfer_matrix{T<:FloatingPoint}(A::Array{T,2},x::Array{T,1},y::Array{T,1},xx::Array{T,1},yy::Array{T,1})

    nr, nc = size(A)

    nx = length(x)
    nxx = length(xx)
    xrange = linrange(x[1],x[end],nx)
    newxrange = linrange(xx[1],xx[end],nxx)

    ny = length(y)
    nyy = length(yy)
    yrange = linrange(y[1],y[end],ny)
    newyrange = linrange(yy[1],yy[end],nyy)

    Ap = zeros(nr,nxx*nyy)

    for i=1:nr
        d = reshape(A[i,:],ny,nx)
        di = Grid.CoordInterpGrid((yrange,xrange),d,Grid.BCnearest,Grid.InterpLinear);
        dp = zeros(nyy,nxx)

        for (xi, xp) in enumerate(newxrange)
            for (yi,yp) in enumerate(newyrange)
                dp[yi,xi] = di[yp,xp]
            end
        end

        Ap[i,:]=reshape(dp,1,int(nxx*nyy))
    end

    return Ap
end

function bilinear(A,x,y,xx,yy)

    nx = length(x)
    xrange = linrange(x[1],x[end],nx)

    ny = length(y)
    yrange = linrange(y[1],y[end],ny)

    Ai = Grid.CoordInterpGrid((yrange,xrange),A,Grid.BCnearest,Grid.InterpLinear)

    nxx = length(xx)
    nyy = length(yy)
    Ap = zeros(nyy,nxx)

    for (xi, xp) in enumerate(xx)
        for (yi,yp) in enumerate(yy)
            Ap[yi,xi] = Ai[yp,xp]
        end
    end

    return Ap
end

function write_to_file(A,energy,pitch,x,sigma,reg_err,dist,dx,dens,dens_err,filename)
    nrows,ncols = size(A)
    npitch,nenergy,nerrs = size(x)

    isfile(filename) && rm(filename)

    nerr_id = NetCDF.NcDim("nerr",nerrs)
    ncols_id = NetCDF.NcDim("ncols",ncols)
    nrows_id = NetCDF.NcDim("nrows",nrows)
    e_id = NetCDF.NcDim("nenergy",nenergy)
    p_id = NetCDF.NcDim("npitch",npitch)

    A_varid = NetCDF.NcVar("transfer_matrix",[nrows_id,ncols_id])
    e_varid = NetCDF.NcVar("energy",e_id)
    p_varid = NetCDF.NcVar("pitch",p_id)
    x_varid = NetCDF.NcVar("x",[p_id,e_id,nerr_id])
    sigma_varid = NetCDF.NcVar("x_err",[p_id,e_id,nerr_id])
    reg_err_varid = NetCDF.NcVar("reg_err",[p_id,e_id,nerr_id])
    dist_varid  = NetCDF.NcVar("dist",[p_id,e_id])
    dx_varid = NetCDF.NcVar("dx",[p_id,e_id,nerr_id])
    dens_varid = NetCDF.NcVar("dens",nerr_id)
    dens_err_varid = NetCDF.NcVar("dens_err",nerr_id)

    ncid =
    NetCDF.create(filename,[A_varid,e_varid,p_varid,x_varid,sigma_varid,reg_err_varid,dist_varid,dx_varid,dens_varid,dens_err_varid], mode = NC_CLASSIC_MODEL)

    NetCDF.putvar(ncid,"transfer_matrix",A)
    NetCDF.putvar(ncid,"energy",energy)
    NetCDF.putvar(ncid,"pitch",pitch)
    NetCDF.putvar(ncid,"x",x)
    NetCDF.putvar(ncid,"x_err",sigma)
    NetCDF.putvar(ncid,"reg_err",reg_err)
    NetCDF.putvar(ncid,"dist",dist)
    NetCDF.putvar(ncid,"dx",dx)
    NetCDF.putvar(ncid,"dens",dens)
    NetCDF.putvar(ncid,"dens_err",dens_err)

    NetCDF.close(ncid)
end
