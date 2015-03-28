

// From the original boost header:

//  stdlibc++ C++0x support is detected via __GNUC__, __GNUC_MINOR__, and possibly
//  __GNUC_PATCHLEVEL__ at the suggestion of Jonathan Wakely, one of the stdlibc++
//  developers. He also commented:
//
//       "I'm not sure how useful __GLIBCXX__ is for your purposes, for instance in
//       GCC 4.2.4 it is set to 20080519 but in GCC 4.3.0 it is set to 20080305.
//       Although 4.3.0 was released earlier than 4.2.4, it has better C++0x support
//       than any release in the 4.2 series."
//
//  Another resource for understanding stdlibc++ features is:
//  http://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html#manual.intro.status.standard.200x
//
// NOTE Mikael Persson:
//  The above is reasonable only under the assumption that GCC is used in conjunction with
//  the libstdc++ library, and with matching versions. What if libstdc++ is used with Clang?

//  C++0x headers in GCC 4.3.0 and later
//
#if !defined( __GLIBCXX__ ) || ( __GLIBCXX__ < 20080305 ) \
  || !( defined( __GXX_EXPERIMENTAL_CXX0X__ ) || ( __cplusplus >= 201103L ) )
#else
#undef BOOST_NO_CXX11_HDR_ARRAY
#undef BOOST_NO_CXX11_HDR_REGEX
#undef BOOST_NO_CXX11_HDR_TUPLE
#undef BOOST_NO_CXX11_HDR_UNORDERED_MAP
#undef BOOST_NO_CXX11_HDR_UNORDERED_SET
#undef BOOST_NO_CXX11_HDR_FUNCTIONAL
#endif

//  C++0x headers in GCC 4.4.0 and later
//
#if !defined( __GLIBCXX__ ) || ( __GLIBCXX__ < 20090421 ) \
  || !( defined( __GXX_EXPERIMENTAL_CXX0X__ ) || ( __cplusplus >= 201103L ) )
#else
#undef BOOST_NO_CXX11_HDR_CONDITION_VARIABLE
#undef BOOST_NO_CXX11_HDR_FORWARD_LIST
#undef BOOST_NO_CXX11_HDR_INITIALIZER_LIST
#undef BOOST_NO_CXX11_HDR_MUTEX
#undef BOOST_NO_CXX11_HDR_RATIO
#undef BOOST_NO_CXX11_HDR_SYSTEM_ERROR
#undef BOOST_NO_CXX11_SMART_PTR
#endif

//  C++0x features in GCC 4.5.0 and later
//
#if !defined( __GLIBCXX__ ) || ( __GLIBCXX__ < 20100414 ) \
  || !( defined( __GXX_EXPERIMENTAL_CXX0X__ ) || ( __cplusplus >= 201103L ) )
#else
#undef BOOST_NO_CXX11_NUMERIC_LIMITS
#undef BOOST_NO_CXX11_HDR_FUTURE
#undef BOOST_NO_CXX11_HDR_RANDOM
#endif

//  C++0x features in GCC 4.6.0 and later
//
#if !defined( __GLIBCXX__ ) || ( __GLIBCXX__ < 20110325 ) \
  || !( defined( __GXX_EXPERIMENTAL_CXX0X__ ) || ( __cplusplus >= 201103L ) )
#else
#undef BOOST_NO_CXX11_HDR_TYPEINDEX
#endif

//  C++0x features in GCC 4.7.0 and later
//
#if !defined( __GLIBCXX__ ) || ( __GLIBCXX__ < 20120322 ) \
  || !( defined( __GXX_EXPERIMENTAL_CXX0X__ ) || ( __cplusplus >= 201103L ) )
#else
#undef BOOST_NO_CXX11_HDR_CHRONO
#undef BOOST_NO_CXX11_ALLOCATOR
#endif

//  C++0x headers not yet (fully!) implemented
//
// GCC 4.8.0 could be __GLIBCXX__ == 20130510
#if !defined( __GLIBCXX__ ) || ( __GLIBCXX__ < 20130510 ) \
  || !( defined( __GXX_EXPERIMENTAL_CXX0X__ ) || ( __cplusplus >= 201103L ) )
#else
#undef BOOST_NO_CXX11_HDR_THREAD
#undef BOOST_NO_CXX11_HDR_TYPE_TRAITS
#undef BOOST_NO_CXX11_HDR_CODECVT
#undef BOOST_NO_CXX11_ATOMIC_SMART_PTR
#endif


// GCC 4.9.0 is __GLIBCXX__ == 20130704


//  --- end ---
