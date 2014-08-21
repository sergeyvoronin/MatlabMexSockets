% wrapper for SVD with socket server; note port is chosen randomly from 1000 to 2000
% please adjust accordingly for different systems
function [U,S,Vt] = svd_wrapper_intel_mkl(A)
    fprintf('picking port..\n');
    port_ok = -1;
    while port_ok < 0
        portno = round(1000 + (2000-1000)*rand);
        fprintf('trying with portno = %d\n', portno);
        command_str = ['lsof -i :', num2str(portno)];
        [status,output] = system(command_str);
        fprintf('status = %d and output = %s\n', status, output);
        if length(output) < 1
            port_ok = 1;
        else
            pause(1.5);
        end
    end

    fprintf('calling socket server..\n');
    command_str = ['./svd_with_socket_server_intel_mkl ', num2str(portno), ' 2>&1 &'];
    fprintf('command_str: %s\n', command_str);
    [status,output] = system(command_str);
    pause(0.75);
    fprintf('status = %d and output = %s\n', status, output);

    fprintf('calling mex_client..\n');
    [U,S,Vt] = svd_mex_client(A,portno);  

    fprintf('ok now kill svd server..\n');
    cmd = ['pkill -9 svd_with_socket_server'];
    system(cmd);
end

