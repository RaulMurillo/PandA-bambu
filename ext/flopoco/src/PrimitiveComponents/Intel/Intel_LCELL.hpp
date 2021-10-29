#ifndef Intel_LCELL_H
#define Intel_LCELL_H

#include "../Primitive.hpp"

namespace flopoco {
///
/// \brief The Intel_LCELL class. Implementation of the target specific Altera lcell_comb.
///
class Intel_LCELL : public Primitive {
        bool _hasSharedArith,_hasDontTouch,_hasLUT7;
      public:

        Intel_LCELL(Operator *parentOp, Target *target, const std::string & lut_mask, const bool &shared_arith = false, const bool &dont_touch = false );

        ~Intel_LCELL() {}

        const bool &hasSharedArith() const{ return _hasSharedArith; }
        const bool &hasLUT7() const{ return _hasLUT7;}
        const bool &hasDontTouch() const{ return _hasDontTouch;}

    };
}//namespace

#endif
