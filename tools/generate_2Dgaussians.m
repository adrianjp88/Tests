function [data] = generate_2Dgaussians(parameters, count, size, noise_type, noise_std_dev)

    fprintf('generating data... ');

    a = parameters.a;
    x0 = parameters.x0;
    y0 = parameters.y0;
    s = parameters.s;
    b = parameters.b;

    if numel(a) < count
        a = repmat(a, 1, count);
    end
    if numel(x0) < count
        x0 = repmat(x0, 1, count);
    end
    if numel(y0) < count
        y0 = repmat(y0, 1, count);
    end
    if numel(s) < count
        s = repmat(s, 1, count);
    end
    if numel(b) < count
        b = repmat(b, 1, count);
    end

    [xi,yj] = meshgrid(0:size-1,0:size-1);

    xi = repmat(xi, [1 1 1])';
    yj = repmat(yj, [1 1 1])';

    data = zeros(size, size, count, 'single');

    for i = 1:count

        if strcmp(noise_type,'gauss') 

            tmp_noise = noise_std_dev * randn(size,size);

            data(:,:,i)...
                = a(i)*exp(-1/2*((xi-x0(i))/s(i)).^2)...
                .*exp(-1/2*((yj-y0(i))/s(i)).^2)...
                + b(i) + tmp_noise;

        elseif strcmp(noise_type,'poisson')

            data(:,:,i)...
                = a(i)*exp(-1/2*((xi-x0(i))/s(i)).^2)...
                .*exp(-1/2*((yj-y0(i))/s(i)).^2)...
                + b(i);

            data(:,:,i) = poissrnd(data(:,:,i),size,size);

        else

            data(:,:,i)...
                = a(i)*exp(-1/2*((xi-x0(i))/s(i)).^2)...
                .*exp(-1/2*((yj-y0(i))/s(i)).^2)...
                + b(i);

        end 

    end

    data = reshape(data,size*size,count);
    
    fprintf('done!\n');

end