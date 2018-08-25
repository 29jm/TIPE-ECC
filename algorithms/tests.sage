# -*- coding: utf-8 -*-
# makes accents work in graphs

# Benchmarks and graphs

load('attacks.sage')
load('cm_approach.sage')
load('rd_approach.sage')
load('schoof.sage')

def test_dlp(dlp, r0, k0, n=10):
    """Solves n ECDLPs using the specified `dlp` algorithm
    at the specified difficulty and returns the average solve time."""
    times = []
    for i in range(n):
        E, r, k, P = gen_random(r0, k0)
        Q = randrange(1, P.order())*P
        t0 = clock()
        m = dlp(P, Q)
        times.append(clock()-t0)
    return sum(times)/n

def macro_test_dlp(dlp, name='DLP solver', n=10):
    """Runs tests for different ECDLP difficulties and creates a graph
    SolveTime = f(difficulty)"""
    bits = []
    times = []
    k0 = 4
    for p_bits in range(10, 26):
        print 'testing with', p_bits, 'bits'
        t = test_dlp(dlp, 2**p_bits, k0, n)
        bits.append(p_bits)
        times.append(t)
        print '->', t, 'seconds on average'
    plt.grid(True)
    plt.xlabel('log_2(|E(F_p)|)')
    plt.ylabel('average time of the '+name)
    plt.title(name+' performance (cofactor %i)' % k0)
    plt.savefig(name+'.png')
    plt.plot(bits, times)
    plt.show()

@parallel(ncpus=4)
def timer(r0, k0, twist):
    from time import clock
    t0 = clock()
    C = gen_random(r0, k0, twist)
    return clock()-t0

def extract(res): # Use with @parallel functions
    return map(lambda x: x[1], res)

def graph_twist_perf(r0, k0):
    N = 100
    perfs = []
    for twist in [False, True]:
        print 'twist:', twist
        times = extract(timer([(r0, k0, twist) for i in range(N)]))
        average = sum(times)/len(times)
        perfs.append(average)

    print 'twists'
    print ' - without:', perfs[0]
    print ' - with:', perfs[1], ((perfs[0]-perfs[1])/perfs[0])*100, '% saved'

    plt.title("Impact du twist sur le temps de generation")
    plt.xlabel("Twist (sans, avec)")
    plt.ylabel("Temps moyen de generation d'une courbe sécurisee (ms)")
    plt.annotate("Curve group size (average): $2^{%i}$" % log(r0, 2),
        xycoords='axes pixels', xy=(10, 340))
    plt.annotate("Curve max cofactor: $%i$" % k0,
        xycoords='axes pixels', xy=(10, 325))
    name = "figures/twist_"+data_name(r0, k0)+".png"
    plt.hist(perfs)
    plt.savefig(name, bbox_inches='tight')
    save_data(r0, k0, perfs)
    plt.clf()

def benchmark(r=2^50, pseudo=False):
    from time import clock

    proof.arithmetic(not pseudo)
    N = 100
    perfs = []

    for i in range(N):
        t0 = clock()
        E, q, k, G = gen_cm(r, 4)
        perfs.append(clock()-t0)

    return perfs

def compare_perfs(r=2^50):
    t1 = mean(benchmark(False, r))
    t2 = mean(benchmark(True, r))

    print 'with proof:', t1
    print 'without proof:', t2
    print 'improvement:', 100*(t1-t2)/t1, '%'

def test_schoof(n=10, r=50):
    """Runs schoof n times with a random curve defined over F_p,
    with p close to r. Returns the timings."""
    from time import clock
    load('rd_approach.sage')
    timings = []
    p = next_prime(r)
    for i in range(n):
        a, b = get_parameters(p)
        t0 = clock()
        t = schoof(a, b, p)
        timings.append(clock()-t0)
        print ' ->', timings[-1]
    return timings

def bench(bits, prefix, gen, args):
    """Saves average compute time for curves with `b` bits of security,
    using the `gen` function to generate curves.
    Results are in 'data/${prefix}_${b}bits_N100.sobj'."""
    import sys
    k0 = 4
    N = 100
    print prefix
    for security in bits:
        print '\nNumber of bits:', security
        r0 = 2^security
        perfs_bits = []
        for i in range(N):
            print i, ' ',
            sys.stdout.flush()
            t0 = clock()
            E, r, k, G = gen(r0, k0, *args)
            perfs_bits.append(clock()-t0)
        name = 'data/'+prefix+'_'+str(security)+'bits_N'+str(N)
        save(perfs_bits, name)

def bench_all():
    bits = range(3, 24+1)
    bench(bits, 'customrd', gen_random, [True])

    bits = range(110, 160, 10)
    bench(bits, 'rd', gen_random)

    bits = range(10, 110, 10) + [140, 160, 180, 224, 256]
    bench(bits, 'cm', gen_cm, [20])

def load_data(bits, prefix, N=100):
    means = []
    devs = []

    for security in bits:
        name = 'data/'+prefix+'_'+str(security)+'bits_N'+str(N)
        data = load(name)

        means.append(mean(data))
        devs.append(std(data))

    return means, devs

def plot_cm():
    import matplotlib.pyplot as plt

    bits = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 140, 180, 224, 256]
    ticks = [10, 30, 50, 70, 100, 140, 180, 224, 256]
    means, devs = load_data(bits, 'cm', 100)

    plt.title(u'Performance de la méthode de la multiplication complexe')
    plt.ylabel(u'Temps moyen de génération (secondes)')
    plt.xlabel(u'Securité de la courbe (bits)')

    plt.xticks(ticks)
    plt.grid(True)
    plt.tight_layout()

    plt.plot(bits, means, 'o-')
    plt.show()

def plot_rd():
    import matplotlib.pyplot as plt

    bits = range(10, 100, 5) + range(100, 150+1, 10)
    ticks = bits
    means, devs = load_data(bits, 'rd', 100)

    plt.title(u'Génération par comptage de points (SEA)')
    plt.ylabel(u'Temps moyen de génération (secondes)')
    plt.xlabel(u'Sécurité de la courbe (bits)')

    plt.xticks(ticks)
    plt.grid(True)
    plt.tight_layout()

    plt.plot(bits, means, 'o-')
    plt.show()

def plot_custom_rd():
    import matplotlib.pyplot as plt

    bits = range(3, 22)
    ticks = bits
    means, devs = load_data(bits, 'customrd', 100)

    plt.title(u'Génération par comptage de points (Schoof)')
    plt.ylabel(u'Temps moyen de génération (secondes)')
    plt.xlabel(u'Sécurité de la courbe (bits)')

    plt.xticks(ticks)
    plt.grid(True)
    plt.tight_layout()

    plt.plot(bits, means, 'o-')
    plt.show()

def plot_both():
    import matplotlib.pyplot as plt

    bits_cm = range(10, 100, 10) + [140, 180]
    bits_rd = range(10, 100, 5) + range(100, 150+1, 10)
    ticks = range(0, 200+1, 20)
    means_rd, devs_rd = load_data(bits_rd, 'rd', 100)
    means_cm, devs_cm = load_data(bits_cm, 'cm', 100)

    plt.title(u'Comparaison des méthodes en performance')
    plt.ylabel(u'Temps moyen de génération (secondes)')
    plt.xlabel(u'Sécurité de la courbe (bits)')

    plt.xticks(ticks)
    plt.grid(True)
    plt.tight_layout()

    plt.plot(bits_rd, means_rd, 'o-')
    plt.plot(bits_cm, means_cm, 'o-')
    plt.show()

def plot_finite_curve():
    E = EllipticCurve(GF(59), [2, 3])
    points = [P.xy() for P in E.points() if P != E(0, 1, 0)]

    x_list = [P[0] for P in points]
    y_list = [P[1] for P in points]

    P, Q = E(20, 45), E(38, 9)
    p, q = P.xy(), Q.xy()
    r = (P+Q).xy()
    s = (-(P+Q)).xy()

    plt.xticks(range(0, 60+1, 10))
    plt.yticks(range(0, 60+1, 10))
    plt.tight_layout()
    plt.title(u"$y^2 = x^3 + 2x + 3$, $p = 59$")

    plt.scatter(x_list, y_list)
    plt.annotate('A', xy=p, xytext=(p[0]+1, p[1]))
    plt.annotate('B', xy=q, xytext=(q[0]+1, q[1]))
    plt.annotate('C', xy=r, xytext=(r[0]+1, r[1]))
    plt.annotate('A+B', xy=s, xytext=(s[0]+1, s[1]))

    plt.show()
