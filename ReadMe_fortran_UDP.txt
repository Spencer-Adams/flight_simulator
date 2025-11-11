The connection type allows for sending/receiving array of values from/to your program. Connections can be made through udp or to a file. The connection type will automatically handle the I/O required for udp or file read/write. The connenction type has an internal timer and will automatically handle the send/receive at the desired rate, so that send/recv commands can be called in the main loop of a program without having to worry about timing.

After initializing, the connection object has either the send command or recv command available, depending on how it was setup. The test.f90 file shows an example of both types.

There are two different udp modules. One for Windows and the other for MacOS. When compiling, just reference the appropriate file. For example, MacOS would use:
gfortran -fdefault-real-8 json_m.f90 jsonx_m.f90 database_m.f90 udp_m.f90 connection_m.f90 test.f90 -o test
while Windows users would use:
gfortran -fdefault-real-8 json_m.f90 jsonx_m.f90 database_m.f90 udp_windows_m.f90 connection_m.f90 test.f90 -o test.exe

Connections are initialized using a JSON file with the following options:

type
    string that is either 'send' or 'receive'
    required

refresh_rate
    refresh rate of the connection in Hz (cyl/sec)
    optional, defaults to 1e12

number_of_values
    integer of the number of values to send/receive
    required

channel_type
    string of the type of connection. It can be 'udp' or 'file'
    required

================ if channel_type is file

filename
    string of the filename
    required

pathname
    string of the path to the file
    optional, will default to the local folder

================ if channel_type is file and type is receive

database_type
    string of the database type. Currently this has to be 'rectilinear'
    required

================ if channel_type is udp

port_ID
    integer of the port number
    required

double_precision
    boolean, if true, values are sent in dp, otherwise they are sent in single precision
    optional, defaults to false

================ if channel_type is udp and type is send

IP_address
    string of the ip address to send to
    optional, defaults to the local machine

================ if channel_type is udp and type is receive

wait_for_input
    boolean, if true, the connection will block code execution and wait for an input whenever recv is called
    optional, defaults to false



