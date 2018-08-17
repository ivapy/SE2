#================================== SE2 ==================================
#
#  function g = SE2(d, rot)
#
#
#  Generates an instance of the class object SE2.  As a class, it acts as
#  though it were a Matlab variable.  Different operations and actions
#  can be applied to it.
#
#
#  Inputs:
#    d		- the displacement/translation.
#    rot	- either the angle of rotation, or the rotation matrix.
#
#================================== SE2 ==================================

classdef(mstring('SE2'), mstring('<handle'))

properties(Constant)
epsilon = 0.0001
end


M()
end

#
methods()

@mfunction("g")
def SE2(d=None, rot=None):

    if (nargin == 0):
        d = mcat([0,0])
        rot = 0
    elif ((size(d, 1) == 1) and (size(d, 2) == 2)):
        d = transpose(d)
    elif ((size(d, 1) != 2) or (size(d, 2) != 1)):
        error(mstring('The translation vector has incorrect dimensions'))
        end

        if (isscalar(rot)):
            rot = mcat([cos(rot), -sin(rot), sin(rot), cos(rot)])
        elif ((size(rot, 1) != 2) or (size(rot, 2) != 2)):
            error(mstring('The rotation input has incorrect dimensions'))
            end

            g.M = mcat([rot, d, 0, 0, 1])

            end

#============================== Adjoint ==============================
 #
 #  function g = adjoint(g1, g2)
 #
 #
 #  Computes and returns the adjoint of g.  When applied to a Lie group
 #  elements, the adjoint is defined to operate as:
 #
 #    Ad_g1 (g2) = g1 * g2 * inverse g1
 #
 #  When applied to a Lie algebra elements in homogeneous form, it
 #  operates as:
 #
 #    Ad_g1 xi_hat = g1 * xi_hat * inverse(g1)
 #
 #  and when applied to a Lie algebra element in vector form it is the
 #  linear matrix operation:
 #
 #                 [ R  | JJ d ]    [ xi1 ]
 #    Ad_g1 xi =   [ --------- ] *  [ xi2 ]
 #                 [ 0  |  1   ]    [ xi3 ]
 #
 #
#============================== Adjoint ==============================

@mfunction("z")
def adjoint(g=None, x=None):

    if (isa(x, mstring('SE2'))):
        z = g * x * inv(g)
    elif ((size(x, 1) == 3) and (size(x, 2) == 1)):
        JJ = mcat([0, 1, OMPCSEMI, -1, 0])
        z = mcat([g.M(mslice[1:2], mslice[1:2]), JJ * g.M(mslice[1:2], 3), OMPCSEMI, 0, 0, 1]) * x
    elif ((size(x, 1) == 3) and (size(x, 2) == 3)):
        z = g.M * x * inv(g.M)
        end

        end

#================================ display ================================
#
#  function display(g)
#
#
#  This is the default display function for the SE2 class.  It simply
#  displays the position followed by the rotation.
#
#================================ display ================================

@mfunction("")
def display(g=None):

    if isequal(get(0, mstring('FormatSpacing')), mstring('compact')):
        disp(mcat([inputname(1), mstring(' =')]))
        disp(g.M)
    else:
        disp(mstring(' '))
        disp(mcat([inputname(1), mstring(' =')]))
        disp(mstring(' '))
        disp(g.M)
        end
        end

#================================ getAngle ===============================
  #
  #  Return the orientation of the Lie group object as an angle.
  #
  #

@mfunction("angle")
def getAngle(g=None):
    angle = atan2(g.M(2, 1), g.M(1, 1))
    end

#============================= getTranslation ============================
  #
  #  Returns the translation or position associated to the group
  #  element.
  #
@mfunction("d")
def getTranslation(g=None):
    d = g.M(mslice[1:2], 3)
    end

#============================== getRotation ==============================
  #
  #  Return the orientation of the Lie group object as a rotation matrix.
  #

@mfunction("R")
def getRotation(g=None):
    R = g.M(mslice[1:2], mslice[1:2])
    end

#================================ inv ================================
  #
  #  invg = inv(g)
  #
  #
  #  Computes and returns the inverse to g.
  #
  #================================ inv ================================

@mfunction("invg")
def inv(g=None):
    tR = transpose(g.M(mslice[1:2], mslice[1:2]))
    invgM = mcat([tR, -tR * g.M(mslice[1:2], 3), OMPCSEMI, 0, 0, 1])
    invg = SE2(invgM(mslice[1:2], end), invgM(mslice[1:2], mslice[1:2]))
    end

#=============================== leftact ==============================
  #
  #  p2 = leftact(g, p)
  #		with p a 2x1 specifying point coordinates.
  #
  #  p2 = leftact(g, v)
  #		with v a 3x1 specifying a velocity.
  #		This applies to pure translational velocities in homogeneous
  #		form, or to SE2 velocities in vector forn.
  #
  #  This function takes a change of coordinates and a point/velocity,
  #  and returns the transformation of that point/velocity under the change
  #  of coordinates.
  #
  #  Alternatively, one can think of the change of coordinates as a
  #  transformation of the point to somewhere else, %  e.g., a displacement
  #  of the point.  It all depends on one's  perspective of the
  #  operation/situation.
  #
  #=============================== leftact ===============================

@mfunction("x2")
def leftact(g=None, x=None):
    if ((size(x, 1) == 2)):
        x2 = g.M * mcat([x, OMPCSEMI, ones(mcat([1, size(x, 2)]))])
        x2 = x2(mslice[1:2], mslice[:])
    elif ((size(x, 1) == 3)):
        L = eye(3)
        L(mslice[1:2], mslice[1:2]).lvalue = g.M(mslice[1:2], mslice[1:2])
        x2 = L * x
        end
        end

#================================== log ==================================
  #
  #  function xi = log(g, tau)
  #
  #  Take the logarithm of the group element g.  If the time period of
  #  the action is not given, it is assumed to be unity.
  #
  #================================== log ==================================

@mfunction("xi")
def log(g=None, tau=None):

    if ((nargin < 2) or isempty(tau)):
        tau = 1
        end

        xi = zeros(mcat([3, 1]))    # Specify size/dimensions of xi.

        #--(1) Obtain the angular velocity.

        xi(3).lvalue = atan2(R(2, 1), R(1, 1)) / tau

        #--(2) Compute the linear velocity.
        if (xi(3) == 0):        # If no rotation, pure translation.
            tau
        else:        # else, use logarithm equations.
            JJ = mcat([0, 1, OMPCSEMI, -1, 0])

            end
            end

#================================ mtimes ===============================
  #
  #  function g = mtimes(g1, g2)
  #
  #
  #  Computes and returns the product of g1 with g2.
  #
  #  Can also be typed as:  >> g3 = g1*g2
  #
  #================================ mtimes ===============================

@mfunction("g")
def mtimes(g1=None, g2=None):
    g = SE2()
    g.M = g1.M * g2.M
    end

#================================== plot =================================
  #
  #  function plot(g, label, linecolor, sc)
  #
  #  Plots the coordinate frame associated to g.  The figure is cleared,
  #  so this will clear any existing graphic in the figure.  To plot on
  #  top of an existing figure, set hold to on.  The label is the name
  #  of label given to the frame (if given is it writen out).  The
  #  linecolor is a valid plot linespec character.  Finally sc is the
  #  specification of the scale for plotting.  It will rescale the
  #  line segments associated with the frame axes and also with the location
  #  of the label, if there is a label.
  #
  #  Inputs:
  #    g		- The SE2 coordinate frame to plot.
  #    label	- The label to assign the frame.
  #    linecolor  - The line color to use for plotting.  (See `help plot`)
  #    sc		- scale to plot things at.
  #		  a 2x1 vector, first element is length of axes.
  #		    second element is a scalar indicating roughly how far
  #		    from the origin the label should be placed.
  #
  #  Output:
  #    The coordinate frame, and possibly a label,  is plotted.
  #
  #================================== plot =================================

@mfunction("")
def plot(g=None, flabel=None, lcol=None, sc=None):

    if ((nargin < 2)):
        flabel = mstring('')
        end

        if ((nargin < 3) or isempty(lcol)):
            lcol = mstring('b')
            end

            if ((nargin < 4) or isempty(sc)):
                sc = mcat([1, 1])
            elif (size(sc, 2) == 1):
                sc = mcat([sc, 2])
                end

                            # Origin of the frame at g's object.

                mcat([sc(1), OMPCSEMI, 0])            # unit axes in object's frame.
                mcat([0, OMPCSEMI, sc(1)])

                isheld = ishold

                lspec = mcat([lcol, mstring('-.')])            # Line specification for plotting.
                pts = mcat([o - 0.5 * ex, o + ex])            # Points of x-axis.
                plot(pts(1, mslice[:]), pts(2, mslice[:]), lspec)
                hold(mstring('on'))
                pts = mcat([o - 0.5 * ey, o + ey])            # Points of y-axis.
                plot(pts(1, mslice[:]), pts(2, mslice[:]), lspec)
                plot(o(1), o(2), mcat([lcol, mstring('o')]), mstring('MarkerSize'), 7)
                # Plot marker at origin.

                if (not isempty(flabel)):                # If label desired, then offset and plot.
                    pts = o - (sc(2) / sc(1)) * (ex + ey) / 4
                    text(pts(1), pts(2), flabel)
                    end

                    if (not isheld):
                        hold(mstring('off'))
                        end

                        axis(mstring('equal'))
                        end

#============================ velocityPlot ===========================
  #
  #  Plots a vector velocity of SE(2) as a vector and a rotation.
  #  Assumes that this is not given in body coordinates, but in the
  #  world frame (so it is in mixed frames for velocities, so to speak).
  #

@mfunction("")
def velocityPlot(g=None, vect=None, rad=None):

    if (nargin < 3):
        rad = 0.5
        end

            # Get base point of twist.

            # Get rotation angle for velocity.
        thetaVals = rotAngle + linspace(0, 3 * pi / 6, 20)
        if (vect(3) < 0):
            thetaVals = -thetaVals
            end
            arcPts = rad * mcat([cos(thetaVals), OMPCSEMI, sin(thetaVals)]) + repmat(basePt, mcat([1, length(thetaVals)]))
            arcVec = diff(arcPts(mslice[:], mslice[end - 1:end]), 1, 2)

            wasHeld = ishold

            hold(mstring('on'))
            if ((vect(1) != 0) or (vect(2) != 0)):
                quiver(basePt(1), basePt(2), vect(1), vect(2), mstring('LineWidth'), 2, mstring('Color'), mcat([0, 0, 1]))
                end

                if (vect(3) != 0):
                    plot(arcPts(1, mslice[1:end - 1]), arcPts(2, mslice[1:end - 1]), mstring('LineWidth'), 2)
                    qh = quiver(arcPts(1, end - 1), arcPts(2, end - 1), arcVec(1), arcVec(2), 2, mstring('LineWidth'), 2, mstring('MaxHeadSize'), 100, mstring('Color'), mcat([0, 0, 1]))
                    end

                    if (not wasHeld):
                        hold(mstring('off'))
                        end

                        end

#============================= twistPlot =============================
  #
  #  Plots a vector velocity of SE(2) as a vector and a rotation at
  #  the object frame defined by g, presuming that the group element
  #  is given in terms of the frame of reference to plot in.  Typically
  #  will be the body frame in the world frame if it is the body
  #  velocity.  Otherwise, should be the identity element to plot as the
  #  spatial velocity.
  #
  #  When rotation is positive, arc will be in first quadrant.
  #  When negative, arc will be in second quadrant.
  #  When zero, no arc.
  #
  #  Likewise, if linear velocity is zero, then no vector.
  #
  #  TODO: Modify so that linespec can be adjusted.
  #

@mfunction("")
def twistPlot(g=None, xi=None, rad=None):

    if (nargin < 3):
        rad = 0.5
        end

            # Get base point of twist.
            # Get rotation matrix to transform ...
        vect = rotMat * xi(mslice[1:2])    #   twist vector to obs frame.

            # Get rotation angle for velocity.
        thetaVals = rotAngle + linspace(0, 3 * pi / 6, 20)
        if (xi(3) < 0):
            thetaVals = -thetaVals
            end
            arcPts = rad * mcat([cos(thetaVals), OMPCSEMI, sin(thetaVals)]); print arcPts
            repmat(basePt, mcat([1, length(thetaVals)]))
            arcVec = diff(arcPts(mslice[:], mslice[end - 1:end]), 1, 2)

            wasHeld = ishold

            hold(mstring('on'))
            if ((xi(1) != 0) or (xi(2) != 0)):
                quiver(basePt(1), basePt(2), vect(1), vect(2), mstring('LineWidth'), 2, mstring('Color'), mcat([0, 0, 1]))
                end

                if (xi(3) != 0):
                    plot(arcPts(1, mslice[1:end - 1]), arcPts(2, mslice[1:end - 1]), mstring('LineWidth'), 2)
                    qh = quiver(arcPts(1, end - 1), arcPts(2, end - 1), arcVec(1), arcVec(2), 0, mstring('LineWidth'), 2, mstring('MaxHeadSize'), 20, mstring('Color'), mcat([0, 0, 1]))
                    get(qh)
                    end

                    if (not wasHeld):
                        hold(mstring('off'))
                        end

                        end

# pole #

@mfunction("qpole")
def pole(g=None):

    end

#=============================== times ===============================
  #
  #  function p2 = times(g, p)
  #
  #
  #  This function is the operator overload that implements the left action
  #  of g on the point p.
  #
  #  Can also be typed as:  >> p2 = g.*p
  #
  #================================= times =================================

@mfunction("p2")
def times(g=None, p=None):
    "p2" = leftact(g, p)
    end

#============================== toHomog ==============================
  #
  #  Return homogeneous form.
  #

@mfunction("M")
def toHomog(g=None):
    "M" = g.M
    end

#============================== toVector =============================
  #
  #  Return vector form.
  #

@mfunction("gvec")
def toVector(g=None):
    gvec = mcat([g.M(mslice[1:2], 3), OMPCSEMI, getAngle(g)])
    end

    end


    methods(Static)

#================================ hat ===============================
  #
  #  Takes a vector form of se(2) and hats it to get the homogeneous
  #  matrix form.
  #

@mfunction("xiHat")
def hat(xiVec=None):
    JJ = mcat([0, 1, OMPCSEMI, -1, 0])
    xiHat = mcat([-JJ * xiVec(3), xiVec(mslice[1:2]), OMPCSEMI, 0, 0, 0])
    end

#=============================== unhat ==============================
  #
  #  Takes a vector form of se(2) and hats it to get the homogeneous
  #  matrix form.
  #

@mfunction("xiVec")
def unhat(xiHat=None):
    xiVec = mcat([xiHat(mslice[1:2], 3), OMPCSEMI, xiHat(2, 1)])
    end

#============================ vec2ProdMat ===========================
  #
  #  Converts the vector form of SE(2) into a producct matrix for use
  #  against vector elements (of points or twists/se(2)).
  #

@mfunction("pMat")
def vec2ProdMat(gVec=None):
    pMat = mcat([cos(g(3)), -sin(g(3)), 0, OMPCSEMI, sin(g(3)), cos(g(3)), 0, OMPCSEMI, 0, 0, 1])
    end

#========================== vec2AdjointMat ==========================
  #
  #  Converts the vector form of SE(2) into an Adjoint matrix for use
  #  against twists/se(2).
  #

@mfunction("AdjMat")
def vec2AdjointMat(gVec=None):
    tVec = mcat([0, 1, OMPCSEMI, -1, 0]) * gVec(mslice[1:2])
    AdjMat = mcat([cos(g(3)), -sin(g(3)), tVec(1), OMPCSEMI, sin(g(3)), cos(g(3)), tVec(2), OMPCSEMI, 0, 0, 1])
    end

# ======================exp==========================
# Computes the exponential of a twist in se(2)
#



@mfunction("gexp")
def exp(xi=None, tau=None):

    if (size(xi, 2) == 1):
        xi = SE2.hat(xi)
        end

        if (nargin < 2):
            tau = 1
            end

            expMat = expm(xi * tau)
            gexp = SE2(expMat(mslice[1:2], 3), expMat(mslice[1:2], mslice[1:2]))

            end


            end

            end








