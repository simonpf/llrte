#ifndef _LLRTE_SURFACES_
#define _LLRTE_SURFACES_

#include "llrte/tracers.h"

namespace llrte::surfaces {

    template<typename Vector>
    class Plane {
    protected:

        Plane(const Vector &base,
              const Vector &normal)
            : base_(base), normal_(normal) {
            // Nothing to do here.
        }

        bool has_crossed(const Vector &v) {
            Vector dv = v - base_;
            return dot(base_, normal_) <= 0.0;
        }

    private:
        Vector base_;
        Vector normal_;
    };

    template <typename Vector>
    class BlackPlane : public Plane<Vector> {

    public:
    using Float = typename Vector::Float;
        BlackPlane(const Vector &base,
                   const Vector &normal)
            : Plane<Vector>(base, normal)
        {
            // Nothing to do here.
        }


    template<typename Photon>
    void apply(Photon &p) {
        absorbed_energy_ += p.get_energy();
        p.set_energy(0.0);
    }


    private:
    Float absorbed_energy_ = 0.0;
    };

} // namespace llrte::surfaces
#endif
