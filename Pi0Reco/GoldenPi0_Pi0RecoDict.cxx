// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GoldenPi0_Pi0RecoDict

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
#include "Pi0Selection.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLPi0Selection_Dictionary();
   static void larlitecLcLPi0Selection_TClassManip(TClass*);
   static void *new_larlitecLcLPi0Selection(void *p = 0);
   static void *newArray_larlitecLcLPi0Selection(Long_t size, void *p);
   static void delete_larlitecLcLPi0Selection(void *p);
   static void deleteArray_larlitecLcLPi0Selection(void *p);
   static void destruct_larlitecLcLPi0Selection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::Pi0Selection*)
   {
      ::larlite::Pi0Selection *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::Pi0Selection));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::Pi0Selection", "Pi0Selection.h", 36,
                  typeid(::larlite::Pi0Selection), DefineBehavior(ptr, ptr),
                  &larlitecLcLPi0Selection_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::Pi0Selection) );
      instance.SetNew(&new_larlitecLcLPi0Selection);
      instance.SetNewArray(&newArray_larlitecLcLPi0Selection);
      instance.SetDelete(&delete_larlitecLcLPi0Selection);
      instance.SetDeleteArray(&deleteArray_larlitecLcLPi0Selection);
      instance.SetDestructor(&destruct_larlitecLcLPi0Selection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::Pi0Selection*)
   {
      return GenerateInitInstanceLocal((::larlite::Pi0Selection*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::Pi0Selection*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLPi0Selection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::Pi0Selection*)0x0)->GetClass();
      larlitecLcLPi0Selection_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLPi0Selection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLPi0Selection(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::Pi0Selection : new ::larlite::Pi0Selection;
   }
   static void *newArray_larlitecLcLPi0Selection(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::Pi0Selection[nElements] : new ::larlite::Pi0Selection[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLPi0Selection(void *p) {
      delete ((::larlite::Pi0Selection*)p);
   }
   static void deleteArray_larlitecLcLPi0Selection(void *p) {
      delete [] ((::larlite::Pi0Selection*)p);
   }
   static void destruct_larlitecLcLPi0Selection(void *p) {
      typedef ::larlite::Pi0Selection current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::Pi0Selection

namespace ROOT {
   static TClass *vectorlElarlitecLcLmcshowergR_Dictionary();
   static void vectorlElarlitecLcLmcshowergR_TClassManip(TClass*);
   static void *new_vectorlElarlitecLcLmcshowergR(void *p = 0);
   static void *newArray_vectorlElarlitecLcLmcshowergR(Long_t size, void *p);
   static void delete_vectorlElarlitecLcLmcshowergR(void *p);
   static void deleteArray_vectorlElarlitecLcLmcshowergR(void *p);
   static void destruct_vectorlElarlitecLcLmcshowergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<larlite::mcshower>*)
   {
      vector<larlite::mcshower> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<larlite::mcshower>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<larlite::mcshower>", -2, "vector", 214,
                  typeid(vector<larlite::mcshower>), DefineBehavior(ptr, ptr),
                  &vectorlElarlitecLcLmcshowergR_Dictionary, isa_proxy, 0,
                  sizeof(vector<larlite::mcshower>) );
      instance.SetNew(&new_vectorlElarlitecLcLmcshowergR);
      instance.SetNewArray(&newArray_vectorlElarlitecLcLmcshowergR);
      instance.SetDelete(&delete_vectorlElarlitecLcLmcshowergR);
      instance.SetDeleteArray(&deleteArray_vectorlElarlitecLcLmcshowergR);
      instance.SetDestructor(&destruct_vectorlElarlitecLcLmcshowergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<larlite::mcshower> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<larlite::mcshower>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElarlitecLcLmcshowergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<larlite::mcshower>*)0x0)->GetClass();
      vectorlElarlitecLcLmcshowergR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElarlitecLcLmcshowergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElarlitecLcLmcshowergR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<larlite::mcshower> : new vector<larlite::mcshower>;
   }
   static void *newArray_vectorlElarlitecLcLmcshowergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<larlite::mcshower>[nElements] : new vector<larlite::mcshower>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElarlitecLcLmcshowergR(void *p) {
      delete ((vector<larlite::mcshower>*)p);
   }
   static void deleteArray_vectorlElarlitecLcLmcshowergR(void *p) {
      delete [] ((vector<larlite::mcshower>*)p);
   }
   static void destruct_vectorlElarlitecLcLmcshowergR(void *p) {
      typedef vector<larlite::mcshower> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<larlite::mcshower>

namespace {
  void TriggerDictionaryInitialization_libGoldenPi0_Pi0Reco_Impl() {
    static const char* headers[] = {
"Pi0Selection.h",
0
    };
    static const char* includePaths[] = {
"/home/david/uboone/larlite/core",
"/home/david/uboone/larlite/UserDev/BasicTool",
"/home/david/SOFTWARE/ROOT_v060410/include",
"/home/david/uboone/larlite/UserDev/GoldenPi0/Pi0Reco/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$Pi0Selection.h")))  Pi0Selection;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "Pi0Selection.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::Pi0Selection", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGoldenPi0_Pi0Reco",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGoldenPi0_Pi0Reco_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGoldenPi0_Pi0Reco_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGoldenPi0_Pi0Reco() {
  TriggerDictionaryInitialization_libGoldenPi0_Pi0Reco_Impl();
}
