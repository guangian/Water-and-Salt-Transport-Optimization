import type { FastifyInstance } from 'fastify';
import type { AppContext } from '../context.js';

export async function healthRoutes(app: FastifyInstance, ctx: AppContext): Promise<void> {
  app.get<{ Querystring: { force?: string } }>('/health', async (req) => ctx.health.report(req.query.force === '1'));
}
