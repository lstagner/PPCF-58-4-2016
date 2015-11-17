function gradEP(E,p; scheme = :Backward, Z = 2.01410178)

    E_dim = length(E)
    p_dim = length(p)

    amu2kg = 1.660538921e-27 # [kg/amu]
    e0 = 1.60217733e-19 # [C]

    m = Z*amu2kg # amu -> kg
    E *= (1e3)*e0 # keV -> joules

    dE = abs(E[2]-E[1])
    dP = abs(p[2]-p[1])

    E_mat = reshape(ones(p_dim)*E',E_dim*p_dim)*ones(E_dim*p_dim)'
    p_mat = reshape(p*ones(E_dim)',E_dim*p_dim)*ones(E_dim*p_dim)'
    
    if scheme == :Backward
        delta_E = zeros(E_dim*p_dim,E_dim*p_dim)
        delta_p = zeros(E_dim*p_dim,E_dim*p_dim)

        for i=1:(E_dim*p_dim)
            if mod(i-1,p_dim) != 0
                delta_p[i,(i-1):i] = [-1.0 1.0]
            end
            if i > p_dim
                delta_E[i,(i-p_dim):(i)] = [-1.0 zeros(1,p_dim-1) 1.0]
            end
        end

    # elseif scheme == :Central
    #     delta_E = diagm([-2*ones(p_dim), zeros((E_dim-2)*p_dim), 2*ones(p_dim)]) .+
    #               diagm([2*ones(p_dim),ones((E_dim-2)*p_dim)],p_dim) .-
    #               diagm([ones((E_dim-2)*p_dim), 2*ones(p_dim)],-p_dim)
    #
    #     delta_E = delta_E/(2*dE)*sqrt(2*m);
    #
    #
    #     Dp_c0 = zeros(E_dim*p_dim)
    #     Dp_c1 = zeros(E_dim*p_dim)
    #     Dp_c2 = zeros(E_dim*p_dim)
    #
    #     for i = 1:E_dim
    #         delta_E[:,(i-1)*p_dim+1:i*p_dim] = delta_E[:,(i-1)*p_dim+1:i*p_dim]*sqrt(E[i])
    #
    #         Dp_c0[(i-1)*p_dim+1:i*p_dim] = [-2 zeros(1,p_dim-2) 2]
    #         Dp_c1[(i-1)*p_dim+1:i*p_dim] = [2 ones(1,p_dim-2) 0]
    #         Dp_c2[(i-1)*p_dim+1:i*p_dim] = [-ones(1,p_dim-2) -2 0]
    #     end
    #
    #     Dp_c1 = Dp_c1[1:end-1]
    #     Dp_c2 = Dp_c2[1:end-1]
    #
    #     delta_p = diagm(Dp_c0) + diagm(Dp_c1,1) + diagm(Dp_c2,-1);
    #     delta_p = delta_p/(2*dP)*sqrt(m/2);
    else
        error("Scheme not implemented")
    end

    J_E = (sqrt(2.0*m)/dE).*sqrt(E_mat)
    J_P = (sqrt(m/2.0)/dP).*sqrt((1.0 .- p_mat.^2)./E_mat)
    J_P[E_mat .== 0.0] = 0.0
    delta_E = J_E.*delta_E
    delta_p = J_P.*delta_p

    return delta_E, delta_p
end
