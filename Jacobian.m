function J = Jacobian(s,d1,d2)
    % construct the jacobian mapping [s,d1,d2] <--> [s,d1+d2,(d1-d2)/s]
    % for the vector of parameter values.
    % s : [N,1]
    % d1: [N,1]
    % d2: [N,1]
    N = numel(s);
    assert((N==numel(d1))&(N==numel(d2)))

    x1 = s;
    x3 = (d2-d1)./s;
    J = zeros([3,3,size(s,3:5)]);
    J(1,1,:) = 2;
    J(2,1,:) = -x3;
    J(3,1,:) = +x3;
    J(2,3,:) = -x1;
    J(3,3,:) = +x1;
    J(2:3,2,:) = 1;
    J = .5*J;
end