import ctypes
lib = ctypes.cdll.LoadLibrary('./librs.so')

F = lib.field_alloc()
gen = lib.polynomial_alloc()

lib.rs_init_field(F, 0x11d)
lib.rs_generator(F, 4, gen)
lib.polynomial_print(gen, "%d")
