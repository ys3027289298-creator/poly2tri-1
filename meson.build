project('poly2tri', ['cpp'])

include = include_directories('.')
lib = static_library('poly2tri', sources : [
	'poly2tri/common/shapes.cc',
	'poly2tri/sweep/advancing_front.cc',
	'poly2tri/sweep/cdt.cc',
	'poly2tri/sweep/sweep.cc',
	'poly2tri/sweep/sweep_context.cc',
])

thread_dep = dependency('threads')
boost_test_dep = dependency('boost', modules : [ 'filesystem', 'unit_test_framework' ], required : false)
if boost_test_dep.found()
	test('Unit Test', executable('unittest', [
		'unittest/main.cpp',
		'unittest/TriangleTest.cpp',
		'unittest/TaubinSmoothingTest.cpp',
	], dependencies : [boost_test_dep, thread_dep],
	include_directories : include_directories('.'),
	cpp_args : ['-DMESH_HAVE_BOOST_TEST'], link_with : lib))
else
	# No Boost: run the self-contained Taubin test through its shim.
	test('Taubin Test', executable('taubin_unittest', [
		'unittest/TaubinSmoothingTest.cpp',
	], include_directories : include_directories('.')))
endif

poly2tri_dep = declare_dependency(include_directories : include, link_with : lib)
