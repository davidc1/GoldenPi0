// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GoldenPi0_TruthFiltersDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "CCpi0Filter.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLCCpi0Filter_Dictionary();
   static void larlitecLcLCCpi0Filter_TClassManip(TClass*);
   static void *new_larlitecLcLCCpi0Filter(void *p = 0);
   static void *newArray_larlitecLcLCCpi0Filter(Long_t size, void *p);
   static void delete_larlitecLcLCCpi0Filter(void *p);
   static void deleteArray_larlitecLcLCCpi0Filter(void *p);
   static void destruct_larlitecLcLCCpi0Filter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::CCpi0Filter*)
   {
      ::larlite::CCpi0Filter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::CCpi0Filter));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::CCpi0Filter", "CCpi0Filter.h", 25,
                  typeid(::larlite::CCpi0Filter), DefineBehavior(ptr, ptr),
                  &larlitecLcLCCpi0Filter_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::CCpi0Filter) );
      instance.SetNew(&new_larlitecLcLCCpi0Filter);
      instance.SetNewArray(&newArray_larlitecLcLCCpi0Filter);
      instance.SetDelete(&delete_larlitecLcLCCpi0Filter);
      instance.SetDeleteArray(&deleteArray_larlitecLcLCCpi0Filter);
      instance.SetDestructor(&destruct_larlitecLcLCCpi0Filter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::CCpi0Filter*)
   {
      return GenerateInitInstanceLocal((::larlite::CCpi0Filter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::CCpi0Filter*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLCCpi0Filter_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::CCpi0Filter*)0x0)->GetClass();
      larlitecLcLCCpi0Filter_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLCCpi0Filter_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLCCpi0Filter(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::CCpi0Filter : new ::larlite::CCpi0Filter;
   }
   static void *newArray_larlitecLcLCCpi0Filter(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::CCpi0Filter[nElements] : new ::larlite::CCpi0Filter[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLCCpi0Filter(void *p) {
      delete ((::larlite::CCpi0Filter*)p);
   }
   static void deleteArray_larlitecLcLCCpi0Filter(void *p) {
      delete [] ((::larlite::CCpi0Filter*)p);
   }
   static void destruct_larlitecLcLCCpi0Filter(void *p) {
      typedef ::larlite::CCpi0Filter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::CCpi0Filter

namespace {
  void TriggerDictionaryInitialization_libGoldenPi0_TruthFilters_Impl() {
    static const char* headers[] = {
"CCpi0Filter.h",
0
    };
    static const char* includePaths[] = {
"/home/david/uboone/larlite/core",
"/home/david/SOFTWARE/ROOT_v60410/include",
"/home/david/uboone/larlite/UserDev/GoldenPi0/TruthFilters/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$CCpi0Filter.h")))  CCpi0Filter;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "CCpi0Filter.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::CCpi0Filter", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGoldenPi0_TruthFilters",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGoldenPi0_TruthFilters_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGoldenPi0_TruthFilters_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGoldenPi0_TruthFilters() {
  TriggerDictionaryInitialization_libGoldenPi0_TruthFilters_Impl();
}
