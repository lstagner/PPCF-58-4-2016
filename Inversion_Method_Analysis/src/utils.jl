using Grid

function resizeTransferMatrix(A,x,y,nxx,nyy)

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
        di = CoordInterpGrid((yrange,xrange),d,BCnearest,InterpLinear);
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

function bilinear(A,x,y,xx,yy)

    nx = length(x)
    xrange = linrange(x[1],x[end],nx)

    ny = length(x)
    yrange = linrange(y[1],y[end],ny)

    Ai = CoordInterpGrid((yrange,xrange),A,BCnearest,InterpLinear)

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
