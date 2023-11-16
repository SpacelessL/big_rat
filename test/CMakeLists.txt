add_definitions(-DBOOST_UT_DISABLE_MODULE)

add_executable(bench_big_int bench_big_int.cpp ../src/big_int.h)
add_executable(ut_big_int ut_big_int.cpp ../src/big_int.h)
add_executable(bench_big_rat bench_big_rat.cpp ../src/big_int.h ../src/big_rat.h)
add_executable(ut_big_rat ut_big_rat.cpp ../src/big_int.h ../src/big_rat.h)
