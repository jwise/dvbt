if bit then
	BOR = bit.bor
	BXOR = bit.bxor
	BAND = bit.band
	BLSH = bit.lshift
	BRSH = bit.rshift
else
	BOR = load("return function(a,b) return a | b end")()
	BXOR = load("return function(a,b) return a ^ b end")()
	BAND = load("return function(a,b) return a & b end")()
	BLSH = load("return function(a,b) return a << b end")()
	BRSH = load("return function(a,b) return a >> b end")()
end

states = {}
-- Generate state table.
for i=0,(BLSH(1, 6) - 1) do
	states[i] = {}
	states[i].val = i
end
for i=0,(BLSH(1, 6) - 1) do
	for inp=0,1 do
		states[i][inp] = {}
		states[i][inp].next = BOR(BRSH(i, 1), BLSH(inp, 5))
		states[i][inp].x = (
			BAND(BRSH(i, 0), 1) +
			BAND(BRSH(i, 3), 1) +
			BAND(BRSH(i, 4), 1) +
			BAND(BRSH(i, 5), 1) +
			inp) % 2
		states[i][inp].y = (
			BAND(BRSH(i, 0), 1) +
			BAND(BRSH(i, 1), 1) +
			BAND(BRSH(i, 3), 1) +
			BAND(BRSH(i, 4), 1) +
			inp) % 2
	end
end

-- Generate bits
fbytes = coroutine.create(function ()
	local f = io.open(arg[1], "rb")
	while true do
		local c = f:read(1)
		if not c then return end
		coroutine.yield(c)
	end
end)

fbits = coroutine.create(function ()
	while true do
		local _,c = coroutine.resume(fbytes)
		if not c then return end
		c = c:byte(1)
		coroutine.yield(BAND(BRSH(c, 7), 1))
		coroutine.yield(BAND(BRSH(c, 6), 1))
		coroutine.yield(BAND(BRSH(c, 5), 1))
		coroutine.yield(BAND(BRSH(c, 4), 1))
		coroutine.yield(BAND(BRSH(c, 3), 1))
		coroutine.yield(BAND(BRSH(c, 2), 1))
		coroutine.yield(BAND(BRSH(c, 1), 1))
		coroutine.yield(BAND(BRSH(c, 0), 1))
	end
end)

-- Initialize the path at step 0.
curst = {}
for sk,_ in pairs(states) do
	curst[sk] = {}
	curst[sk].st = sk
	curst[sk].prev = nil
	curst[sk].pm = 0
end

-- Inductive step functiom from n -> n + 1.
function viterbi(curst, x, y)
	local newst = {}
	for sk,_ in pairs(states) do
		for inp=0,1 do
			-- compute a branch metric
			local bm = 0
			if x and states[sk][inp].x ~= x then bm = bm + 1 end
			if y and states[sk][inp].y ~= y then bm = bm + 1 end
			
			-- compute a path metric
			local pm = curst[sk].pm + bm
			
			-- if it's the shortest path metric for the state we land in, update
			local next = states[sk][inp].next
			if not newst[next] or newst[next].pm > pm then
				newst[next] = {}
				newst[next].st = sk
				newst[next].pm = pm
				newst[next].prev = { st = curst[sk], inp = inp }
			end
		end
	end
	return newst
end

while true do
	local alive,x = coroutine.resume(fbits)
	if not alive then print(x) abort() end
	if not x then break end

	local alive,y = coroutine.resume(fbits)
	if not alive then print(y) abort() end
	if not y then break end
	
	curst = viterbi(curst, x, y)
	
	local alive,y = coroutine.resume(fbits)
	if not alive then print(y) abort() end
	if not y then break end
	
	curst = viterbi(curst, nil, y)
end

-- Now do the reverse pass.
local best = nil
for sk,_ in pairs(states) do
	if not best or best.pm < curst[sk].pm then
		best = curst[sk]
	end
end

bits = {}

io.stderr:write(string.format("Final path metric was %d.\n", best.pm))
while best and best.prev do
	-- io.stderr:write(string.format("pm %d, inp %d, st %x\n", best.pm, best.prev.inp, best.st))
	table.insert(bits, best.prev.inp)
	best = best.prev.st
end
io.stderr:write(string.format("Initial state was %x.\n", best.st))

for i=1,#bits,8 do
	if #bits - i < 8 then break end
	local c = 0
	c = BOR(c, BLSH(bits[#bits - i + 1 - 0], 7))
	c = BOR(c, BLSH(bits[#bits - i + 1 - 1], 6))
	c = BOR(c, BLSH(bits[#bits - i + 1 - 2], 5))
	c = BOR(c, BLSH(bits[#bits - i + 1 - 3], 4))
	c = BOR(c, BLSH(bits[#bits - i + 1 - 4], 3))
	c = BOR(c, BLSH(bits[#bits - i + 1 - 5], 2))
	c = BOR(c, BLSH(bits[#bits - i + 1 - 6], 1))
	c = BOR(c, BLSH(bits[#bits - i + 1 - 7], 0))
	io.stdout:write(string.char(c))
end
