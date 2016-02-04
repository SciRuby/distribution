# Added by John O. Woods, SciRuby project.
module Distribution
  module Gamma
    module Ruby_
      class << self
        include Math

        # Gamma random number generator
        # == Arguments
        # * +a+: alpha or k
        # + +b+: theta or 1/beta
        #
        # Adapted the function itself from GSL-2.1 in rng/gamma.c: gsl_ran_gamma
        #
        # ==References
        # * http://www.gnu.org/software/gsl/manual/html_node/The-Gamma-Distribution.html
        # * http://en.wikipedia.org/wiki/Gamma_distribution
        def rand(a,b)
          if a < 1
            u = Distribution::Normal.rng.call
            x = rand(1.0 + a, b) * u ** (1.0 / a)
            return x
          else
            d = a - 1.0 / 3.0
            c =  (1.0 / 3.0) / Math.sqrt(d)
            v = 0
            while 1

              while v <= 0
                x = Distribution::Normal.rng.call
                v = 1.0 + c * x
              end

              v = v * v * v
              u = Distribution::Uniform.rng.call
              if u < 1 - 0.0331 * x * x * x * x
                break
              end
              if log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))
                break
              end
            end
            return b * d * v
          end
        end

        # Return a Proc object which returns a random number drawn
        # from the Gamma random number generator for given (alpha,beta)
        def rng(a,b)
          -> { rand(a, b)}
        end

        # Gamma distribution probability density function
        #
        # If you're looking at Wikipedia's Gamma distribution page, the arguments for this pdf function correspond
        # as follows:
        #
        # * +x+: same
        # * +a+: alpha or k
        # + +b+: theta or 1/beta
        #
        # This is confusing! But we're trying to most closely mirror the GSL function for the gamma distribution
        # (see references).
        #
        # Adapted the function itself from GSL-1.9 in rng/gamma.c: gsl_ran_gamma_pdf
        #
        # ==References
        # * http://www.gnu.org/software/gsl/manual/html_node/The-Gamma-Distribution.html
        # * http://en.wikipedia.org/wiki/Gamma_distribution
        def pdf(x, a, b)
          return 0 if x < 0
          if x == 0
            return 1.quo(b) if a == 1
            return 0
          elsif a == 1
            Math.exp(-x.quo(b)).quo(b)
          else
            Math.exp((a - 1) * Math.log(x.quo(b)) - x.quo(b) - Math.lgamma(a).first).quo(b)
          end
        end

        # Gamma cumulative distribution function
        def cdf(x, a, b)
          return 0.0 if x <= 0.0

          y = x.quo(b)
          return (1 - Math::IncompleteGamma.q(a, y)) if y > a
          (Math::IncompleteGamma.p(a, y))
        end
      end
    end
  end
end
