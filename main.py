from PIL import Image

img1 = Image.open('Images/old.jpg~')
img2 = Image.open('Images/new.jpg~')

ls1 = [((114, 154), (152, 153)), ((192, 157), (233, 159)), ((104, 135),
(130, 130)), ((131, 131), (154, 137)), ((194, 139), (213, 135)), ((215, 135),
(238, 138)), ((170, 143), (173, 197)), ((149, 208), (169, 215)),
((172, 217), (194, 208)), ((135, 238), (170, 240)), ((171, 240), (205, 238)),
((148, 169), (119, 220)), ((195, 170), (221, 215)), ((177, 3), (126, 11)),
((126, 13), (95, 49)), ((94, 49), (76, 90)), ((74, 94), (87, 161)), ((87, 163),
(92, 213)), ((92, 214), (106, 261)), ((108, 263), (161, 289)), ((178, 2),
(247, 27)), ((247, 29), (273, 73)), ((273, 74), (278, 116)), ((278, 116),
(261, 168)), ((260, 167), (256, 219)), ((255, 220), (245, 255)), ((245, 255),
(184, 288)), ((106, 263), (107, 296)), ((108, 298), (171, 358)), ((245, 257),
(246, 294)), ((246, 297), (189, 359)), ((78, 144), (72, 176)), ((72, 175),
(87, 213)), ((280, 150), (285, 183)), ((285, 183), (263, 219))]

ls2 = [((127, 141), (168, 144)), ((204, 142), (239, 145)), ((121, 125),
(147, 118)), ((149, 118), (174, 131)), ((194, 129), (215, 122)),((218, 122),
(245, 127)), ((186, 133), (190, 187)), ((162, 193), (186, 203)), ((190, 202),
(205, 192)), ((152, 217), (186, 224)), ((188, 225), (222, 219)), ((170, 155),
(136, 197)), ((204, 155), (250, 190)), ((194, 9), (130, 23)), ((128, 24),
(101, 70)), ((101, 72), (97, 104)), ((97, 107), (109, 155)), ((109, 158),
(113, 206)), ((113, 208), (130, 256)), ((132, 259), (179, 274)), ((197, 10),
(262, 30)), ((263, 32), (287, 71)), ((286, 74), (292, 107)), ((292, 108),
(273, 161)), ((273, 163), (270, 208)), ((270, 209), (261, 244)), ((260, 247),
(203, 275)), ((130, 264), (128, 292)), ((128, 296), (199, 366)), ((261, 249),
(267, 289)), ((268, 293), (207, 366)), ((102, 137), (95, 172)), ((97, 175),
(110, 204)), ((292, 145), (297, 173)), ((297, 176), (273, 209))]

neighbor = lambda m, x, y: m[int(x), int(y)]

def bilinear(m, x, y):
	x1 = int(x)
	x2 = x1 + 1
	y1 = int(y)
	y2 = y1 + 1

	wx1, wx2 = x - x1, x2 - x
	wy1, wy2 = y - y1, y2 - y

	try:
		a = zip(map(lambda c: wx1 * c, m[x1, y1]), map(lambda c: wx2 * c, m[x2, y1]))
		a = map(lambda t: t[0] + t[1], a)

		b = zip(map(lambda c: wx1 * c, m[x1, y2]), map(lambda c: wx2 * c, m[x2, y2]))
		b = map(lambda t: t[0] + t[1], b)

		r = zip(map(lambda c: wy1 * c, a), map(lambda c: wy2 * c, b))
		return map(lambda t: int(t[0] + t[1]), r)

		# a = wx1 * m[x1, y1] + wx2 * m[x2, y1]
		# b = wx1 * m[x1, y2] + wx2 * m[x2, y2]
		# return int(wy1 * a + wy2 * b + 0.5)

	except (IndexError, ValueError):
		return neighbor(m, x, y)

def reverse_map(original, new_size, inverse_transform, resample):
	new = Image.new(original.mode, new_size, (255, 255, 255))
	new_pixels = new.load()
	original_pixels = original.load()
	
	max_x, max_y = new_size

	for x in range(max_x):
		for y in range(max_y):
			u, v = inverse_transform(x, y)

			try:
				new_pixels[x, y] = resample(original_pixels, u, v)

			except (IndexError, ValueError):
				pass

	return new

def cross_dissolve(original1, original2, alpha):
	max_x, max_y = original1.size
	result = Image.new(original1.mode, (max_x, max_y))

	result_pixels = result.load()
	original1_pixels = original1.load()
	original2_pixels = original2.load()

	for x in range(max_x):
		for y in range(max_y):
			complement = 1 - alpha
			temp = zip(map(lambda c: complement * c, original1_pixels[x, y]), map(lambda c: alpha * c, original2_pixels[x, y]))
			result_pixels[x, y] = map(lambda t: int(t[0] + t[1]), temp)
			
			# result_pixels[x, y] = int((1 - alpha) * original1_pixels[x, y] + alpha * original2_pixels[x, y] + 0.5)

	return result

def warp(line_pairs):
	from math import sqrt

	def vectorize(p1, p2):
		x1, y1 = p1
		x2, y2 = p2

		return x2 - x1, y2 - y1

	def add(*vs):
		sx = sy = 0

		for x, y in vs:
			sx += x
			sy += y
		
		return sx, sy

	def scale(s, v):
		x, y = v

		return s * x, s * y

	def internal_multiply(v1, v2):
		x1, y1 = v1
		x2, y2 = v2

		return x1 * x2 + y1 * y2

	def perpendicular(v):
		a, b = v

		if a == 0:
			return b, 0

		x = -(b / a)
		y = sqrt( internal_multiply(v, v) / (x * x + 1) )
		x = x * y

		return x, y

	def pre_calculate(line_pair):
		(ps, qs), (pd, qd) = line_pair

		pdqd = vectorize(pd, qd)
		aux = internal_multiply(pdqd, pdqd)
		npdqd = sqrt(aux)
		pdqdp = perpendicular(pdqd)

		psqs = vectorize(ps, qs)
		psqsp = perpendicular(psqs)
		npsqs = sqrt( internal_multiply(psqs, psqs) )

		return pd, qd, ps, pdqd, pdqdp, aux, npdqd, psqs, npsqs, psqsp

	A = 0.1
	B = 1.5
	P = 1
	cs = map(pre_calculate, line_pairs)

	def transform(dx, dy):
		d = (dx, dy)

		def diplacement(c):
			pd, qd, ps, pdqd, pdqdp, aux, npdqd, psqs, npsqs, psqsp = c
			pdd = vectorize(pd, d)
			u = internal_multiply(pdd, pdqd) / aux
			v = internal_multiply(pdd, pdqdp) / npdqd
			s = add(ps, scale(u, psqs), scale(v / npsqs, psqsp))

			if 0 < u < 1:
				dist = abs(v)
			elif u < 0:
				dist = sqrt( internal_multiply(pdd, pdd) )
			else:
				qdd = vectorize(qd, d)
				dist = sqrt( internal_multiply(qdd, qdd) )

			return vectorize(d, s), dist

		ds, dists = zip( *map(diplacement, cs) )
		ws = map(lambda p: ( p[0][8] ** P / (A + p[1]) ) ** B, zip(cs, dists))
		dsum = add( *map(lambda p: scale(p[0], p[1]), zip(ws, ds)) )
		wsum = reduce(lambda x, y: x + y, ws)

		return add(d, scale(1 / wsum, dsum))

	return transform

def main(n):
	def interpolate(a):
		c = 1 - a

		def aux(lp):
			(p1, q1), (p2, q2) = lp
			pi = (c * p1[0] + a * p2[0], c * p1[1] + a * p2[1])
			qi = (c * q1[0] + a * q2[0], c * q1[1] + a * q2[1])

			return pi, qi

		return map(aux, lps)

	lps = zip(ls1, ls2)

	for i in range(1, n):
		a = float(i) / n
		lsi = interpolate(a)
		warped1 = reverse_map(img1, img1.size, warp( zip(ls1, lsi) ), bilinear)
		warped2 = reverse_map(img2, img1.size, warp( zip(ls2, lsi) ), bilinear)
		cross_dissolve(warped1, warped2, a).save( 'Images/' + str(i) + '.jpg' )

if __name__ == '__main__':
	from sys import argv

	if len(argv) > 1: 
		main(int(argv[1]))