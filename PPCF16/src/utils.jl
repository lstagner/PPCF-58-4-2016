using Grid

function resize_transfer_matrix{T<:AbstractFloat}(A::AbstractArray{T,2},
                                                  x::AbstractArray{T,1},
                                                  y::AbstractArray{T,1},
                                                  nxx::Int,nyy::Int)

    nr, nc = size(A)

    nx = length(x)
    xrange = range(x[1],(x[end]-x[1])/nx,nx)
    newxrange = range(x[1],(x[end]-x[1])/nxx,nxx)

    ny = length(y)
    yrange = range(y[1],(y[end]-y[1])/ny,ny)
    newyrange = range(y[1],(y[end]-y[1])/nyy,nyy)

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

        Ap[i,:]=reshape(dp,1,Int(nxx*nyy))
    end

    return Ap,[newxrange],[newyrange]
end

function resize_transfer_matrix{T<:AbstractFloat}(A::AbstractArray{T,2},
                                                  x::AbstractArray{T,1},
                                                  y::AbstractArray{T,1},
                                                  xx::AbstractArray{T,1},
                                                  yy::AbstractArray{T,1})

    nr, nc = size(A)

    nx = length(x)
    nxx = length(xx)
    xrange = range(x[1],(x[end]-x[1])/nx,nx)
    newxrange = range(xx[1],(xx[end]-xx[1])/nxx,nxx)

    ny = length(y)
    nyy = length(yy)
    yrange = range(y[1],(y[end]-y[1])/ny,ny)
    newyrange = range(yy[1],(yy[end]-yy[1])/nyy,nyy)

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

        Ap[i,:]=reshape(dp,1,Int(nxx*nyy))
    end

    return Ap
end

function bilinear(A,x,y,xx,yy)

    nx = length(x)
    xrange = range(x[1],(x[end]-x[1])/nx,nx)

    ny = length(y)
    yrange = range(y[1],(y[end]-y[1])/ny,ny)

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
