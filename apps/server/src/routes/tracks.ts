import type { FastifyInstance } from 'fastify';
import type { AppContext } from '../context.js';
import { ApiError } from '../types.js';

export async function trackRoutes(app: FastifyInstance, ctx: AppContext): Promise<void> {
  app.get('/tracks', async () => ctx.tracks.list());

  app.get('/knowledge', async () => ctx.knowledge.list());

  app.get<{ Params: { slug: string } }>('/knowledge/:slug', async (req) => {
    const slug = decodeURIComponent(req.params.slug ?? '');
    if (!slug) throw new ApiError(400, 'bad_slug', 'slug 不能为空');
    return ctx.knowledge.read(slug);
  });
}
