
gfortran -fdefault-real-8 -O2 json.f90 jsonx.f90 linalg.f90 database_m.f90 upd_windows_m.f90 connection_m.f90 adams.f90 timing.f90 sim.f90 main.f90 -o sim.exe -lws2_32 