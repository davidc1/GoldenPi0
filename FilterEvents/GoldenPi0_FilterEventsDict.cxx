// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GoldenPi0_FilterEventsDict

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
#include "FilterCCpi0.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLFilterCCpi0_Dictionary();
   static void larlitecLcLFilterCCpi0_TClassManip(TClass*);
   static void *new_larlitecLcLFilterCCpi0(void *p = 0);
   static void *newArray_larlitecLcLFilterCCpi0(Long_t size, void *p);
   static void delete_larlitecLcLFilterCCpi0(void *p);
   static void deleteArray_larlitecLcLFilterCCpi0(void *p);
   static void destruct_larlitecLcLFilterCCpi0(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::FilterCCpi0*)
   {
      ::larlite::FilterCCpi0 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::FilterCCpi0));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::FilterCCpi0", "FilterCCpi0.h", 25,
                  typeid(::larlite::FilterCCpi0), DefineBehavior(ptr, ptr),
                  &larlitecLcLFilterCCpi0_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::FilterCCpi0) );
      instance.SetNew(&new_larlitecLcLFilterCCpi0);
      instance.SetNewArray(&newArray_larlitecLcLFilterCCpi0);
      instance.SetDelete(&delete_larlitecLcLFilterCCpi0);
      instance.SetDeleteArray(&deleteArray_larlitecLcLFilterCCpi0);
      instance.SetDestructor(&destruct_larlitecLcLFilterCCpi0);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::FilterCCpi0*)
   {
      return GenerateInitInstanceLocal((::larlite::FilterCCpi0*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::FilterCCpi0*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLFilterCCpi0_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::FilterCCpi0*)0x0)->GetClass();
      larlitecLcLFilterCCpi0_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLFilterCCpi0_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLFilterCCpi0(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::FilterCCpi0 : new ::larlite::FilterCCpi0;
   }
   static void *newArray_larlitecLcLFilterCCpi0(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::FilterCCpi0[nElements] : new ::larlite::FilterCCpi0[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLFilterCCpi0(void *p) {
      delete ((::larlite::FilterCCpi0*)p);
   }
   static void deleteArray_larlitecLcLFilterCCpi0(void *p) {
      delete [] ((::larlite::FilterCCpi0*)p);
   }
   static void destruct_larlitecLcLFilterCCpi0(void *p) {
      typedef ::larlite::FilterCCpi0 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::FilterCCpi0

namespace {
  void TriggerDictionaryInitialization_libGoldenPi0_FilterEvents_Impl() {
    static const char* headers[] = {
"FilterCCpi0.h",
0
    };
    static const char* includePaths[] = {
"/uboone/app/users/davidc1/LArLite/core",
"/grid/fermiapp/products/larsoft/root/v6_04_02/Linux64bit+2.6-2.12-e7-prof/include",
"/uboone/app/users/davidc1/LArLite/UserDev/GoldenPi0/FilterEvents/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$FilterCCpi0.h")))  FilterCCpi0;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "FilterCCpi0.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::FilterCCpi0", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGoldenPi0_FilterEvents",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGoldenPi0_FilterEvents_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGoldenPi0_FilterEvents_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGoldenPi0_FilterEvents() {
  TriggerDictionaryInitialization_libGoldenPi0_FilterEvents_Impl();
}
