function init(obj, symbolic_x, symbolic_f, symbolic_g)
    if isempty(symbolic_x) || isempty(symbolic_f) || isempty(symbolic_g)
        error('x,f,g is empty. Create a class funtcion defineSystem and define your dynamics symbolically.')
    end
    
    if ~isa(symbolic_f, 'sym')
        f_ = sym(symbolic_f);
    else
        f_ = symbolic_f;
    end

    if ~isa(symbolic_g, 'sym')
        g_ = sym(symbolic_g);
    else
        g_ = symbolic_g;
    end

    x           = symbolic_x;
    obj.xdim    = size(x,1);
    obj.udim    = size(g_,2);
    obj.f       = matlabFunction(f_,'vars',{x});
    obj.g       = matlabFunction(g_,'vars',{x});
    
    % Taylor expand the system 
    u           = sym('u',[obj.udim 1]);
    dx          = obj.f(x) + obj.g(x)*u;
    dx_app      = vpa(taylor(dx,x,'Order',obj.paras.r));
    obj.fpoly   = matlabFunction(dx_app, 'vars', {x,u});

    % Linearize system matrices
    obj.A       = matlabFunction(jacobian(dx,x), 'vars', {x});
    obj.B       = matlabFunction(jacobian(dx,u), 'vars', {x});
end