function follow_contour(f, p; boundary = x -> true, step_range = (0.01, 0.1), tol = 0.001, orien = 0, maxiter=5000, verbose=true)
    boundary(p) || error("Starting point outside boundary")

    fc = x -> (f(x) - f(p))^2

    hmin, hmax = step_range

    ri = copy(p)
    curve = Polygon()
    sinth = sin(0.1)
    costh = cos(0.1)

    res = ForwardDiff.HessianResult(p)
    cnt = 0
    while cnt < maxiter
        if !boundary(ri)
            verbose && info("Hit Boundary at ", ri)
            return Nullable{Polygon{Float64}}()
        end

        if fc(ri) > tol
            #println(cnt, " cor ", ri)
            ForwardDiff.hessian!(res,fc,ri)
            grad = ForwardDiff.gradient(res)
            hess = ForwardDiff.hessian(res)
            dr = hess\grad
            ri = ri - dr
            cnt=cnt+1
            continue
        end
        #println(cnt, " ", ri)
        push!(curve.vertices,(ri[1],ri[2]))

        ForwardDiff.hessian!(res,f,ri)
        ci = ForwardDiff.value(res)
        grad = ForwardDiff.gradient(res)
        hess = ForwardDiff.hessian(res)
        gn = norm(grad)
        gn <= 0.0 && error("Zero Gradient encountered")

        g = grad/gn
        t = [-g[2],g[1]]

        k = -(t'*hess*t)[1]/(gn)
        rc = 1/k

        #if hmin < 0.1*abs(rc) < hmax
        #    rp = ri + rc*g - abs(rc)*costh*g + abs(rc)*sinth*t
        #else
        #    h = clamp(0.1*abs(rc),hmin,hmax)
        #    th = min(pi/4,abs(h/rc))
        #    rp = ri + rc*g - abs(rc)*cos(th)*g + abs(rc)*sin(th)*t
        #end

        h = clamp(abs(0.1*rc),hmin,hmax)
        rp = ri + h*t

        dp = norm(p-rp)
        di = norm(p-ri)
        d = norm(rp-ri)

        if ((dp^2 + di^2 - d^2)/(2*di*dp)) < 0.0
            push!(curve.vertices,(p[1],p[2]))
            break
        end

        ri = rp
        cnt = cnt+1
    end

    if length(curve.vertices) == 1
        push!(curve.vertices,(p[1],p[2]))
    end

    if cnt >= maxiter && verbose
        warn("Maximum Iteration exceeded: $maxiter ", p)
    end

    if orien != 0 && orien*orientation(curve) < 0
        reverse!(curve.vertices)
    end
    return Nullable{Polygon{Float64}}(curve)
end
