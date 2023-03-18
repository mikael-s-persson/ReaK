#Keep it Hot !The secret to high - performance code

This tutorial is the first time I write about coding techniques for improving the performance of software. This is going to be a step-by-step journey through a real example, with real benchmark results at every step so that you can see the effect of the different optimization techniques.

And as my first tutorial on performance optimization, I have to talk about the most important principle in performance optimization: keeping things hot. If you expect this tutorial to be about showing bits of code and the assembly listings generated by the compiler and then, showing how to reduce the number of instructions using nifty tricks, then you are sorely mistaken about what has most impact on performance.

When I named this tutorial "Keep it Hot", I wasn't thinking about the [Cameo song from the 80s](https://www.youtube.com/watch?v=uLl4_rNuK6s), although you're welcome to listen to it while you read this tutorial, that is, if you agree with me that a funky beat is the *real* secret to producing good code!

You might be wondering what it is exactly that is supposed to be kept *hot*. Well, software is about using memory and machine instructions to compute values. So, what should be kept hot? Values, memory and machine instructions. And what does it mean to be "hot"? It means to be in high demand, as in the expression "a hot commodity".

## Background on the real example

Let me just tell you how this tutorial came about. I wrote some simple collision detection code a while back and never got around to optimizing its performance because it was good enough as it was and I had greater worries elsewhere. But two things were clear: (1) I knew exactly what I had to do to optimize it, and (2) when I would optimize it, I would record benchmarks at every step and write a tutorial detailing the process I went through. And luckily (or "as expected", if I wasn't so humble), it all worked out and you are reading the result right now.

So, the case study for this tutorial is a piece of collision detection code. It's a simple set of functions to find collisions between a limited set of 3D primitive bounding volumes: spheres, boxes, planes, cylinders and capped-cylinders (or "capsule" or "pill", i.e., a cylinder with hemi-spherical ends). The idea with these primitive shapes is that you can wrap almost any weird object with a set of primitive shapes to create a kind of invisible protective volume around it. This is sometimes used in physics engines and dynamic simulators (although more precise contact points are typically required), but I use it for detecting proximity to obstacles when planning motions to be executed by a robot. In that application domain, obstacle avoidance, there is no need for anything more exact, because being conservative is the name of the game.

With 5 types of shapes, it requires a total of 15 individual functions to compute the minimum distance or contact point for each possible combination of shapes. Most of those functions can be written using some simple math (geometry), the sphere-to-sphere test being the simplest of such *closed-form* calculations. There are a few combinations that are a bit more difficult and require the use of an iterative algorithm, and I just used the GJK-EPA algorithm because it's one of the fastest and it's very easy to write.

This example is also interesting to look at because it's a classic *double-dispatch* problem. If you're not familiar with the term, it means that unlike the usual *single-dispatch* situation where each derived class implements its own virtual functions, you have a situation where the choice of function to be called depends on two classes, not just one. In other words, given two objects, each being of a particular shape, the decision of which function to call depends on both types of shapes.

## The original solution

Obviously, I cannot post the entire code here, because it's far too much code. But like I said at the start, this is not about nitty-gritty coding details, because they don't really matter that much to the performance. So, I'm just going to give you a bird-eyes view of it. The full code can be found in [my repository](https://github.com/mikael-s-persson/ReaK/tree/modular/libs/geometry) (but note that some details are omitted or changed in the code presented in this tutorial, just to keep things cleaner).

Each primitive shape is derived from a base class called `shape_3D`. The shape classes themselves are just called with the very unoriginal names `sphere`, `box`, `plane`, `cylinder` and `capped_cylinder`. They each contain the information necessary to define their dimensions. More importantly, they are all positioned using a position and orientation relative to some "anchor" coordinate system, which would typically be the body-fixed frame of the object they belong to, and often, those frames form a chain of relative transformations all the way to the global coordinate system.

To construct a complete model, many shapes have to be assembled together. The class responsible for this is called `proxy_query_model_3D` (I used the term *proximity query* instead of collision detection, because my code also produces minimum distance calculations, not just collisions), and here are the relevant parts of it:

    class proxy_query_model_3D {
 private:
  std::vector<std::shared_ptr<shape_3D>> mShapeList;

 public:
  friend class proxy_query_pair_3D;

  proxy_query_model_3D& addShape(const std::shared_ptr<shape_3D>& aShape);

  // ...
    };

Easy enough so far, right?

The friend declaration in the above is for a second class that is used when two of these models are involved in a proximity query, and its relevant parts are as follows:

    class proxy_query_pair_3D {
 private:
  std::shared_ptr<proxy_query_model_3D> mModel1;
  std::shared_ptr<proxy_query_model_3D> mModel2;

 public:
  void setModelPair(const std::shared_ptr<proxy_query_model_3D>& aModel1,
                    const std::shared_ptr<proxy_query_model_3D>& aModel2);

  proximity_record_3D findMinimumDistance() const;

  bool gatherCollisionPoints(std::vector<proximity_record_3D> & aOutput) const;

  // ..
    };

The `findMinimumDistance` and `gatherCollisionPoints` are the two functions that do all the interesting stuff, that is, compute the minimum distance between the two models or gather all the collision points between them, respectively. The `proximity_record_3D` class simply contains result information.

The `gatherCollisionPoints` function implementation basically looks like this:

    bool proxy_query_pair_3D::gatherCollisionPoints(std::vector< proximity_record_3D >& aOutput) const {
  for (auto& pshape1 : mModel1->mShapeList) {
    vect<double, 3> p1 = pshape1->getPose().getGlobalPose().Position;
    for (auto& pshape2 : mModel2->mShapeList) {
      vect<double, 3> p2 = pshape2->getPose().getGlobalPose().Position;

      if (norm_2(p2 - p1) - pshape1->getBoundingRadius() -
              pshape2->getBoundingRadius() >
          0.0)
        continue;

      proximity_record_3D tmp = compute_proximity(*pshape1, *pshape2);
      if (tmp.mDistance < 0.0)
        aOutput.push_back(tmp);
    };
  };
  return !aOutput.empty();
    };

The idea there is that the nested loops traverse through all the possible pairs of shapes between the two models. Then, an early check based on their bounding radius and relative distance is made to avoid computing the proximity between them if they cannot possibly collide. And then, the `compute_proximity` function is called, which is the function that performs the double dispatch. And finally, if a collision occurred, signaled by a negative minimum distance between them, then the proximity record is appended to the output vector.

The double-dispatch problem is notorious because there isn't really a nice solution to it. There are a few nifty tricks to clean the syntax (one of which I will almost accidentally fall into later in this tutorial), but the basic pattern to do it is like this:

    proximity_record_3D compute_proximity(const shape_3D& shape1, const shape_3D& shape2) {
  if (typeid(shape1) == typeid(box)) {
    const box& b1 = static_cast<const box&>(shape1);
    if (typeid(shape2) == typeid(box)) {
      const box& b2 = static_cast<const box&>(shape2);
      return compute_proximity(
          b1, b2);  // call the box-to-box overload of 'compute_proximity'
    } else if (typeid(shape2) == typeid(sphere)) {
      const sphere& s2 = static_cast<const sphere&>(shape2);
      return compute_proximity(
          b1, s2);  // call the box-to-sphere overload of 'compute_proximity'
    }               // .. so on..
  } else if (typeid(shape1) == typeid(sphere)) {
    const sphere& s1 = static_cast<const sphere&>(shape1);
    if (typeid(shape2) == typeid(box)) {
      const box& b2 = static_cast<const box&>(shape2);
      return compute_proximity(
          s1, b2);  // call the sphere-to-box overload of 'compute_proximity'
    }               // ... so on...
  }  // .... so on, so forth... I know, this is tedious and ugly!!
    };

There are overloads of the `compute_proximity` function for each possible combination of shapes. Those are the functions that actually compute the proximity / collision points.

This is essentially all there is to it, at least, all you really need to know about it for now.

Finally, to benchmark the performance of this code, I created a simple program that generates two models containing a large number of primitives whose placement and dimensions are all randomly generated, and then, the calls to the `gatherCollisionPoints` function are timed reliably by calling it a large number of times for a large number of randomly generated models. At the end, the benchmark spits out the average time per call to the `gatherCollisionPoints` function. With the original implementation as described above, the performance obtained was as follows:

    N=50     620,000 ns
    N=100  2,550,000 ns

Which means that when each model contains 50 shapes, it takes 620 microseconds to gather all the collision points between them. And when each model contains 100 shapes, it takes 2.5 milliseconds. 

This is our starting point. Ready to optimize?

N.B.: Before someone points out that I seem to be abusing the use of `shared_ptr`, just note that I am aware of this (and it is indeed abusive), but there are reasons for that, and I don't want to get into that here.

## Keep values hot!

The first low hanging fruit for optimization is the pre-computation of values that are in high demand during the proximity computations. In particular, the computation of the global poses of the shapes, that is, the position and orientation of the shape in the global coordinate system, which is computed by the `getGlobalPose` function in the nested loops above. This global pose is a per-shape computation, not a per-pair computation, and therefore, it is wasteful to compute it inside those nested loops. Furthermore, it is very likely (and in fact, known for sure) to be required to be computed inside the individual `compute_proximity` functions for each shape, and that calculation is clearly redundant to the one done in the loops.

Pre-computing those poses can be done quite simply as so:

    bool proxy_query_pair_3D::gatherCollisionPoints(std::vector< proximity_record_3D >& aOutput) const {
  std::vector<pose_3D<double>> mdl1_poses, mdl2_poses;
  for (auto& pshape1 : mModel1->mShapeList)
    mdl1_poses[&pshape1 - mModel1->mShapeList.data()] =
        pshape1->getPose().getGlobalPose();
  for (auto& pshape2 : mModel2->mShapeList)
    mdl2_poses[&pshape2 - mModel2->mShapeList.data()] =
        pshape2->getPose().getGlobalPose();

  for (auto& pshape1 : mModel1->mShapeList) {
    const vect<double, 3>& p1 =
        mdl1_poses[&pshape1 - mModel1->mShapeList.data()].Position;
    for (auto& pshape2 : mModel2->mShapeList) {
      const vect<double, 3>& p2 =
          mdl2_poses[&pshape2 - mModel2->mShapeList.data()].Position;

      if (norm_2(p2 - p1) - pshape1->getBoundingRadius() -
              pshape2->getBoundingRadius() >
          0.0)
        continue;

      proximity_record_3D tmp = compute_proximity(
          *pshape1, mdl1_poses[&pshape1 - mModel1->mShapeList.data()], *pshape2,
          mdl2_poses[&pshape2 - mModel2->mShapeList.data()]);

      if (tmp.mDistance < 0.0)
        aOutput.push_back(tmp);
    };
  };
  return !aOutput.empty();
    };

Where you can notice that the pre-computed poses are handed to the `compute_proximity` function call so that they can be used within them, eliminating the need to recompute them inside those functions.

We could go further with this and see what data could be pre-computed for each kind of shapes, but this will be sufficient for now.

Running the benchmark program again produces the following results:

    N=50     140,000 ns
    N=100    540,000 ns

What a great start! We just went from 620ms to 140ms, and from 2,550ms to 540ms. That's more than 75% faster!

## Keep instructions hot!

Now, we'll get into an issue that might not be obvious to a lot of people. There is a function that I didn't put up yet because it was pretty obvious what it did, this is the `addShape` function from the `proxy_query_model_3D` class. Here is the original definition of it:

    proxy_query_model_3D& proxy_query_model_3D::addShape(const std::shared_ptr< shape_3D >& aShape) {
  mShapeList->push_back(aShape);
  return *this;
    };

Pretty obvious, isn't it? But there is one important optimization we can make here, and the beauty of it is that it will only take one additional line of code:

    proxy_query_model_3D& proxy_query_model_3D::addShape(const std::shared_ptr< shape_3D >& aShape) {
  mShapeList.push_back(aShape);
  std::inplace_merge(mShapeList.begin(), mShapeList.end() - 1, mShapeList.end(),
                     [](const std::shared_ptr<shape_3D>& lhs,
                        const std::shared_ptr<shape_3D>& rhs) {
                       return typeid(*lhs).before(typeid(*rhs));
                     });
  return *this;
    };

Ok... so, I had to spread that line over 4 lines... but hey, it's still one expression. What is happening here is that I sort (incrementally, using the standard `inplace_merge` algorithm) the shapes in the list by their derived type (from `typeid` and the `before` function from the standard `type_info` class). If you are unfamiliar with the `[]() { }` construct that you see above, that is a lambda expression, which is a feature that was introduced in C++11 to create a simple function object on site. Note that this could also be done with a separate function (to be compatible with C++98/03) as so:

    static bool compare_shape_types(const std::shared_ptr< shape_3D >& lhs, const std::shared_ptr< shape_3D >& rhs) {
  return typeid(*lhs).before(typeid(*rhs));
    };
    
    proxy_query_model_3D& proxy_query_model_3D::addShape(const std::shared_ptr< shape_3D >& aShape) {
  mShapeList.push_back(aShape);
  std::inplace_merge(mShapeList.begin(), mShapeList.end() - 1, mShapeList.end(),
                     compare_shape_types);
  return *this;
    };

Running the benchmark program again produces the following results:

    N=50     122,000 ns
    N=100    480,000 ns

This doesn't look as dramatic as the previous step, but we are still looking at about a 10% improvement in performance from the previous benchmark results, which is pretty significant considering that we haven't touched any of the code inside of the loops that compute the proximity. What's going on here?

What's happening here is that you have to understand that all the instructions that make up your program are stored in sections that we call "functions", and these instructions need to exist in memory somewhere, and more importantly, they have to be delivered to the CPU to be executed. When your program runs, the entire program is loaded in RAM memory, but RAM is very very slow relative to the CPU clock. So, to speed things up, memory will be temporarily stored on the CPU cache(s), which is much faster. The section of cache that stores instructions is called the *instruction cache*, or commonly called *I-cache*. On top of that, modern CPUs have deep instruction pipelines that line up many instructions in advance to be able to stream line and possibly interleave their execution. 

What all of this boils down to is that whenever you call a function, one of three things can happen: 

 1. The function call could be so predictable that the instructions of that function are already in the pipeline before you even call it.
 2. The instructions of that function had already been cached before, and therefore, the function call is a simple matter of grabbing the instructions from the CPU cache and line them up to be put in the pipeline, which takes some time, but it's not too bad, especially if you hit L1 cache (it'll be worse if you hit L2 or L3 cache).
 3. The instructions of that function are not in the CPU cache, which means that they have to be fetched all the way back to RAM... and that's slow, real slow, like hundreds or thousands of clock-cycles slow. We call that a *cache miss*.

So, the point is that anything that has been called recently in the past is likely to still be on the CPU cache, i.e., those are "hot functions". And, you want to keep hot functions hot as long as possible and avoid calling "cold" functions. And that's where the sorting of the shapes comes into play. By sorting shapes by type, boxes are with boxes, spheres are with spheres, and so on. When we traverse the sets of shapes during the nested loops in the `gatherCollisionPoints` function, each consecutive call to compute proximity is probably going to involve the same types of shapes, and thus, involve the "hot" functions that were just called on the last computation. This reduces the I-cache misses and improves the performance.

The impact of this was not that great in this case, but believe me, this can bring major performance benefits, especially on tighter loops, that is, when the number of instructions done within the iterations are not as many. In other words, this is a matter of comparative costs between cache misses (100-300 clock cycles) and instructions to execute (nominally, 1 clock cycle per instruction). So, a single cache miss is equivalent to hundreds of wasted clock cycles that you could have spent on more useful things. What our benchmark results tell us is that in our previous version of the code, 10% of the time was spent waiting on "cold" instructions to be transfered from RAM (or L3/L2 cache) to L1 cache.

## Keep memory hot!

What is good for the goose is good for the gander. All these CPU cache issues that I just talked about apply to memory (data)
just as well as they apply to instructions.So,
    your data should also be kept hot to try to streamline the use of the CPU
        cache.

    The main problem here is with the `mShapeList` vector :

    std::vector<std::shared_ptr<shape_3D>>
        mShapeList;

This is a vector of pointers to shapes.The problem with this is that whenever we
    iterate over this vector(as we do a lot),
    what we really want to do is iterate over the shapes,
    but we are doing so *indirectly *.Moreover,
    those shapes are each dynamically allocated individually,
    meaning that they are likely to be scattered all over RAM memory.So,
    every time we go to fetch a shape,
    it is likely to reside in some random piece of
    "cold" memory somewhere in RAM,
    and it will need to be brought,
    at a snail's pace, to the CPU cache where it can be used.

    So,
    what we really should have is a vector like this :

    std::vector<shape_3D> mShapeList;

But there's a problem here. The `shape_3D` is a base-class, and our original `shared_ptr<shape_3D>` pointers are actually pointing to different kinds of shapes. Here is how we can get out of this conundrum, with a nice little class template called `boost::variant` (from the [Boost.Variant](http://www.boost.org/doc/libs/1_57_0/doc/html/variant.html) library):

    typedef boost::variant<box, sphere, plane, cylinder, capped_cylinder>
        any_shape;

std::vector<any_shape> mShapeList;

What the `boost::variant` class does is create a `union` of all the types given
    to it,
    and additionally,
    it also identifies which type of object it currently holds,
    and is thus called a *discriminating union
            *.The consequence of this is that all the shapes are now stored by
                value inside the `mShapeList` vector
                    .This will make the traversals more efficient in terms of
                        the locality of memory access patterns
                    .

        The `boost::variant` class template is used in a somewhat peculiar
            way.It uses a kind of generic visitation pattern,
    which is actually pretty nice,
    but it's not something I want to delve into too much here. Suffice to say, here is how the final `gatherCollisionPoints` function looks:

    namespace {  // anonymous namespace to hide things from the linker.

  struct convert_to_shape_3D_visitor {
    shape_3D& operator()(shape_3D& s) const { return s; };
    typedef shape_3D& result_type;
  };

  struct compute_proximity_3D_visitor {
    const pose_3D<double>* p1;
    const pose_3D<double>* p2;
    compute_proximity_3D_visitor(const pose_3D<double>& aP1,
                                 const pose_3D<double>& aP2)
        : p1(&aP1), p2(&aP2){};
    template <typename T, typename U>
    proximity_record_3D operator()(const T& s1, const U& s2) const {
      return compute_proximity(s1, *p1, s2,
                               *p2);  // double-dispatch by overloading here.
    };
    typedef proximity_record_3D result_type;
  };

  struct get_bounding_radius_3D_visitor {
    double operator()(const shape_3D& s) const {
      return s.getBoundingRadius();
    };
    typedef double result_type;
  };
};

bool proxy_query_pair_3D::gatherCollisionPoints(
    std::vector<proximity_record_3D>& aOutput) const {
  std::vector<pose_3D<double>> mdl1_poses, mdl2_poses;
  for (auto& vshape1 : mModel1->mShapeList)
    mdl1_poses[&vshape1 - mModel1->mShapeList.data()] =
        boost::apply_visitor(convert_to_shape_3D_visitor(), vshape1)
            .getPose()
            .getGlobalPose();
  for (auto& vshape2 : mModel2->mShapeList)
    mdl2_poses[&vshape2 - mModel2->mShapeList.data()] =
        boost::apply_visitor(convert_to_shape_3D_visitor(), vshape2)
            .getPose()
            .getGlobalPose();

  for (auto& vshape1 : mModel1->mShapeList) {
    const vect<double, 3>& p1 =
        mdl1_poses[&vshape1 - mModel1->mShapeList.data()].Position;
    for (auto& vshape2 : mModel2->mShapeList) {
      const vect<double, 3>& p2 =
          mdl2_poses[&vshape2 - mModel2->mShapeList.data()].Position;

      if (norm_2(p2 - p1) -
              boost::apply_visitor(get_bounding_radius_3D_visitor(), vshape1) -
              boost::apply_visitor(get_bounding_radius_3D_visitor(), vshape2) >
          0.0)
        continue;

      proximity_record_3D tmp = boost::apply_visitor(
          compute_proximity_3D_visitor(
              mdl1_poses[&vshape1 - mModel1->mShapeList.data()],
              mdl2_poses[&vshape2 - mModel2->mShapeList.data()]),
          vshape1, vshape2);

      if (tmp.mDistance < 0.0)
        aOutput.push_back(tmp);
    };
  };
  return !aOutput.empty();
};

I agree that this is not pretty, but optimized code is rarely pretty. But I should note that this way of doing things removes the need for the long and tedious double-dispatching `compute_proximity` function (the one with all the nested if-else-if blocks). This is because the multi-visitor version of the `boost::apply_visitor` function of the Boost.Variant library will extract the type of both variant objects and call the visitor with both objects resolved to their actual types, which means that the call to `compute_proximity` that is done inside the `compute_proximity_3D_visitor` will actually get dispatched to the correct proximity computation algorithm for the particular pair of types involved, through overloading. This is pretty much the best way to solve the double-dispatch problem, as far as I know.

Running the benchmark program again produces the following results:

    N=50     106,000 ns
    N=100    405,000 ns

This is about 16% improvement over the previous version, and again, that is without really changing any of the core code that actually does the calculations.

You might be skeptical about this way of doing things, because after all, this seems to completely defeat the purpose of using base classes and polymorphism, which is true. But the purpose of software is not to glorify and worship the almighty object-oriented design patterns, its purpose is to process data. This philosophy is called *data-oriented design* and is a major driving principle in high-performance computing, especially in the computer game's industry.

## Conclusion

It's over already, and we have a final version of the code that is more than 6 times faster than the original one. And I must remind you that this was done without looking at any assembly listings, and thank goodness, without having to write any assembly. It did not require using any special tools like profilers or special optimizers, like whole-program optimizers or profile-guided optimizers (btw, I tried both of those to see, and they had no effect).

It's really important to understand those issues and how to keep things hot, because there is no amount of microscopic code-tinkering that will solve any of these issues, and in general, profilers and special optimizers will not be able to help you either. Generally-speaking, you can trust the compiler to do a good job with microscopic optimizations, and they do have a significant impact on performance, but compilers are really good at optimizing that stuff.

But the kinds of techniques that I demonstrated in this tutorial are not things that any compiler would ever be able to do automatically. You, the programmer, determine the data layout, the (pre-)computation of values, and the sequence of function calls. If you don't streamline them for performance, nobody will, and no tool will do it for you.

And if you take this case study as an example, notice that the individual collision detection functions probably include a few thousands instructions each. Optimizing these instructions by inspecting the assembly listings might allow you to squeeze a hand-full of instructions out of them, but that could take you hours of analysis to figure out and you will only end up with about 1% improvement in performance. In some rare circumstances, this 1% might be worth your time, but you will get a lot more mileage by first applying the principles that I showed in this tutorial.

## Going beyond

The principles discussed here generalize far beyond the scope of the case study used here. For instance, keeping computed values hot is a very important thing to keep an eye on when writing code. There are a lot of situations in which it pays off to do pre-computations or splitting an algorithm into parts that need to belong inside the tight loop and parts that can be factored out as a pre-computation step. There are also general techniques that relate to this principle. For instance, *memoization* is about automatically storing previous results of a certain complex computation to avoid having to re-compute them whenever the computation is invoked with the same input data. There is also a technique for the inverse case, which is called *lazy evaluation*, which is about delaying complex computations as much as possible to avoid computing them until the very moment the results of them are needed. So, the point is that if values are in high-demand, pre-compute them, and if they are in low-demand, compute them at the last minute.

Optimizing for instruction cache misses is a much harder problem in general, and the compiler can do some of this work with profile-guided optimization, although this is very tricky to use and get good results with. At the very least, doing something like what I have demonstrated here, that is, sorting your objects or tasks such that the same code is repeatedly invoked, is a very easy way to get some significant improvements in performance without having to dig too deep into the assembly code or the specifics of your target platform. Also, if you are thinking about using parallelism in your code, these types of techniques will become very important.

And finally, optimizing memory layouts is one of the most important optimizations you can make. In this particular case study, it wasn't as significant as I had hoped, in part because it wasn't a very data-intensive task. In tasks that are more data-intensive, optimizing your memory access patterns is extremely important. In those types of situation, it is not uncommon to get between 10x and 100x performance improvements from this issue alone. Moreover, cache misses get *exponentially worse* as the amount of data grows. This is because as the overall memory used by your program increases, the likelihood that any arbitrary memory address is cold increases, which means that performance degrades very rapidly. For example, it is not uncommon for large binary search trees to exhibit a linear performance curve, instead of the theoretical logarithmic performance curve, purely because of cache misses. Very often, memory access patterns are easy to optimize when you know what to look for, as was the case here, but sometimes it is not as easy, especially with dynamic programming problems and graph algorithms, in which case, you might want to look into *cache-oblivious* algorithms and data-structures which are designed to produce a near-optimal average-case cache performance on *any* CPU architecture.