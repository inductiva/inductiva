## Cost Optimization Strategies

### 💰 Spot Instances: Up to 90% Savings

**How it works:** Use spare cloud capacity at heavily discounted rates
**Risk:** Machines can be preempted (reclaimed) with 30-second notice
**Best for:** Fault-tolerant simulations, jobs <12 hours

**Preemption probability by duration:**
- <2 hours: ~5% chance
- 2-6 hours: ~15% chance  
- 6-12 hours: ~25% chance
- >24 hours: ~60% chance

**Spot strategy tips:**
1. **Checkpoint frequently** to enable easy restarts
2. **Use multiple regions** to reduce preemption risk
3. **Monitor spot history** in our [machine catalog](https://inductiva.ai/machines)
4. **Test resilience** with intentional interruptions

### 🕐 Right-Sizing Your Machines

**Common oversizing mistakes:**
- Using highmem when standard suffices
- Choosing latest generation when older works fine
- Over-provisioning vCPUs for single-threaded code

**Right-sizing checklist:**
- [ ] Memory usage stays below 80% of allocation
- [ ] CPU utilization averages 60-90%
- [ ] Simulation completes without OOM errors
- [ ] No idle cores during computation phases

### 📊 Cost Monitoring

Track your spending patterns:
- **Per-simulation costs** via Inductiva dashboard
- **Resource utilization** through system metrics
- **Spot vs. regular pricing** trends
- **Regional price differences**

## Advanced Selection Strategies

### 🔄 Multi-Machine Workflows

Sometimes multiple smaller machines outperform one large machine:

**When to use multiple machines:**
- Embarrassingly parallel problems
- Independent parameter sweeps
- Fault tolerance requirements
- Cost optimization for long runs

**Example:** Instead of c3d-standard-96 → 3x c3d-standard-32
- **Benefit:** Better fault tolerance, potential cost savings
- **Drawback:** More complex orchestration

### 🌍 Regional Considerations

Machine availability and pricing vary by region:
- **US regions:** Best availability, competitive pricing
- **European regions:** GDPR compliance, slightly higher costs
- **Asian regions:** Lower latency for Asia-based teams
- **Check availability:** Some machine types limited to specific regions

### 🔧 Custom Requirements

**Network-intensive simulations:**
- Choose machines with higher network bandwidth
- Consider placement groups for multi-machine setups

**Storage-intensive workloads:**
- Attach high-performance SSDs
- Consider local SSD vs. persistent disk trade-offs

**Compliance requirements:**
- Confidential computing options available
- Specific regions for data residency